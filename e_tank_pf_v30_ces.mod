%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% E-TANK PERFECT FORESIGHT -- V30-CES (Production finale CES gigogne)
%%%% Auteure : Asmae EL MOUHSSINE -- Juin 2026
%%%%
%%%% STRUCTURE DE PRODUCTION FINALE (CES gigogne, suggestion Schubert) :
%%%%   V_KL  = (K_y)^alpha_KL * L_a^(1-alpha_KL)                    (composite K-L, Cobb-Douglas interne)
%%%%   Y_g   = ass * [ m_Y*V_KL^(1-1/sigma_Y) + (1-m_Y)*E_y^(1-1/sigma_Y) ]^(1/(1-1/sigma_Y))
%%%%   La Cobb-Douglas (benchmark historique) est le cas limite sigma_Y -> 1,
%%%%   avec m_Y = alpha_y+beta_y = 0.93 et alpha_KL = alpha_y/(alpha_y+beta_y) = 0.3226
%%%%   (valide par Test 1 : ecarts < 0.01pp vs CD a sigma_Y=0.9999).
%%%%
%%%% DIFFERENCES V30 vs V29 (4 changements chirurgicaux, herites de la base CD) :
%%%%   (1) Energie propre : E_c = A_c * K_c^alpha_c * L_u_c^(1-alpha_c)
%%%%       alpha_c = 0.20 -- corrige l'inversion de prime salariale
%%%%       P_c = Cobb-Douglas unit cost ; K = K_y + K_c (equilibre capital)
%%%%   (2) sigma_L = 1.40 (vs 0.80) -- Katz-Murphy / Acemoglu-Autor
%%%%   (3) omega   = 0.012 (vs 0.07) -- burden_H ~ 0.19
%%%%   (4) ebar    = 0.044 (vs 0.05) -- plancher subsistance
%%%%
%%%% VARIABLES ENDOGENES : 47 (46 base V30 + V_KL, composite CES gigogne)
%%%% EQUATIONS : 47 (Marche des biens implicite par Walras)
%%%% VAREXO : tau, P_oil, e_ROW
%%%%
%%%% EQUILIBRE GENERALE :
%%%%   - Deux types de menages : HtM (lambda=0.21) + Ricardiens (1-lambda=0.79)
%%%%   - Production finale CES gigogne (K-L)/Energie ; shares au SS calibres pour matcher la CD
%%%%   - Recyclage taxe carbone : R_tau = tau*E_d
%%%%       theta*R_tau -> subvention R&D vert
%%%%       phi*(1-theta)*R_tau -> transferts HtM
%%%%       (1-phi)*(1-theta)*R_tau -> transferts Ricardiens
%%%%   - Fermeture open-economy : exports = facture petrole (tb=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; warning off;

%%
%%%%%%%%%%%%% Variables endogenes (45) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var
% Preferences Stone-Geary CES (3)
P_Z         % indice de prix CES energie/non-energie
Z_H         % composite bien-etre HtM (Stone-Geary CES)
Z_R         % composite bien-etre Ricardien

% Consommation par type et bien (4)
C_H_n       % conso non-energie HtM
C_R_n       % conso non-energie Ricardien
C_H_e       % conso energie HtM (= ebar + supranumeraire)
C_R_e       % conso energie Ricardien

% Offre de travail (2)
L_H_u       % heures travaillees HtM (non-qualifie)
L_R_s       % heures travaillees Ricardien (qualifie)

% Decisions ricardiens : capital, investissement, MU (4)
Lambda_R    % utilite marginale revenu Ricardien
I           % investissement (detrended)
K           % capital (detrended, stock fin de periode)
q_K         % q de Tobin (Ricardien)

% Prix et quantites energie (retailer CES + mix friction) (9)
E           % energie totale (composite CES)
E_c         % energie propre (domestique)
E_d         % energie sale (importee)
E_y         % energie utilisee en production
m           % ratio mix E_c/E_d (variable d'etat : m(-1) dans FOC)
P_e         % prix de l'energie composite (multiplicateur Lagrange retailer)
P_c         % prix energie propre = w_u/A_c
P_d         % prix energie sale = er*P_oil + tau
er          % taux de change reel

% Production (10)
Y_g         % output brut (avant dommages)
Y_net       % output net de dommages = (1-Dam)*Y_g
r_k         % taux de location du capital
L_a         % agregat CES de travail
L_s_y       % travail qualifie en production finale
L_u_y       % travail non-qualifie en production finale
L_u_c       % travail non-qualifie en energie propre (emplois verts)
K_c         % capital dans le secteur energie propre (V30 : K = K_y + K_c)
L_R_c       % travail qualifie en R&D vert
w_s         % salaire qualifie
w_u         % salaire non-qualifie
W_a         % salaire moyen de l'agregat CES
V_KL        % composite capital-travail de la production finale (CES gigogne)

% Innovation & croissance (2)
A_c         % stock de technologie verte (stationnaire, EE : delta_c*A_c_ss = eta_c*L_R_c_ss)
gamma_agg   % taux de croissance semi-endogene (reporting : gz*Y_net/Y_net(-1))

% Climat (3)
e_em        % emissions domestiques = phi_d*E_d
X           % stock de carbone atmospherique
Dam         % dommages climatiques Dam = d0+d1*X+d2*X^2

% Recyclage taxe carbone (4)
R_tau       % recettes taxe carbone = tau*E_d
T_H         % transferts vers HtM
T_R         % transferts vers Ricardiens
C_n         % agregat conso non-energie = lambda*C_H_n + (1-lambda)*C_R_n

% Commerce exterieur (1)
EX          % exportations (paient la facture petroliere)

% Bien-etre (recursions, pour CEV) (2)
W_H         % bien-etre HtM (somme actualisee)
W_R         % bien-etre Ricardien
;

%%
%%%%%%%%%%%%% Variables exogenes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo
tau         % taxe carbone (trajectoire deterministe Net-Zero)
P_oil       % prix mondial du petrole (AR1 ou trajectoire)
e_ROW       % emissions reste du monde (exogene)
;

%%
%%%%%%%%%%%%% Parametres %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters
% Preferences
beta        % facteur d'escompte
sigma       % CRRA (EIS = 1/sigma)
psi         % elasticite Frisch inverse
chi         % poids du travail (calibre : L_ss = 1/3)
lambda_htM  % part HtM (Campbell-Mankiw)
omega       % poids energie dans CES (budget share energie)
epsilon_E   % elasticite substitution energie/non-energie (< 1 : energie rigide)
ebar        % plancher subsistance energie (Stone-Geary)

% Production KLE (parts somment a 1)
alpha_y     % part capital (= 0.30)
beta_y      % part travail (= 0.63)
% part energie = 1-alpha_y-beta_y = 0.07 (implicite)
alpha_c     % part capital dans energie propre (V30 : Cobb-Douglas K_c^alpha_c * L_u_c^(1-alpha_c))
sigma_L     % elasticite substitution qualifie/non-qualifie (Katz-Murphy)
omega_L     % poids qualifie dans CES travail (calibre : w_s/w_u = 1.60)
delta       % taux de depreciation capital
phi_k       % cout d'ajustement investissement (I/K ratio)
ass         % TFP neutrale (calibree pour Y_ss = 1)

% Energie / mix
alpha_E     % poids energie propre dans CES energie
rho_E       % elasticite substitution propre/sale
phi_mix     % friction d'ajustement du mix (greenflation)

% Innovation & croissance
gz          % taux de croissance BGP (frontiere mondiale, = 1.003)
delta_c     % taux de depreciation technologie verte
eta_c       % productivite R&D vert (calibre : delta_c*A_c_ss = eta_c*L_R_c_ss)

% Climat
phi_d       % emissions par unite energie sale
delta_clim  % absorption carbone atmospherique
d0 d1 d2   % parametres fonction dommages (DICE quadratique)

% Open economy
eta_X       % elasticite prix exportations (Henriet et al.)
D_ex_ss     % demande etrangere (calibree pour er_ss = 1)

% Politique
theta       % part R_tau vers R&D (vs transferts) -- levier 1
phi         % progressivite transferts (part HtM) -- levier 2

% Etat stationnaire (pour initval + diagnostics)
gz_ss       % = gz (identite)
tau_ss      % taxe carbone initiale (= tau_init)
P_oil_ss    % prix petrole SS
A_c_ss      % techno verte SS

% CES production finale (suggestion Schubert)
alpha_KL    % part capital dans composite CD interne V = Kprod^alpha_KL * L_a^(1-alpha_KL)
m_Y         % poids du composite KL dans la CES externe (calibre au SS)
sigma_Y     % elasticite substitution (KL)/energie dans la CES externe
;

%%
%%%%%%%%%%%%% SECTION MATLAB : calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% PARAMETRES STRUCTURELS (directement du tableau V7)
% =========================================================================
beta       = 0.997;      % discount factor
sigma      = 1.50;       % CRRA (EIS = 0.67)
psi        = 2.0;        % inverse Frisch
lambda_htM = 0.21;       % part HtM
omega      = 0.012;      % energie weight in CES (V30 : reduit pour burden_H ~ 0.19)
epsilon_E  = 0.30;       % energie/non-energie substitution (< 1)
alpha_y    = 0.30;       % capital share
beta_y     = 0.63;       % labour share
% energy share = 1 - 0.30 - 0.63 = 0.07 (vérification: somme = 1 ✓)
alpha_c    = 0.20;       % V30 : part capital dans energie propre (Cobb-Douglas)
sigma_L    = 1.40;       % V30 : skilled/unskilled substitution (Katz-Murphy : 1.4-2.0)
delta      = 0.025;      % depreciation trimestrielle ~10%/an
phi_k      = 4.0;        % investment adjustment cost
alpha_E    = 0.25;       % clean weight in energy CES
rho_E      = 1.80;       % clean/dirty substitution
phi_mix    = 2.0;        % mix friction (greenflation)
gz         = 1.003;      % BGP growth (~1.2%/yr quarterly)
gz_ss      = gz;
delta_c    = 0.005;      % knowledge obsolescence quarterly (~2%/an)
phi_d      = 1.0;        % emissions per unit dirty energy (normalise)
delta_clim = 0.0035;     % carbon absorption (DICE)

% --- CES production finale (suggestion Schubert) ---
sigma_Y    = 0.50;                        % benchmark : complementarite (KL)/energie
alpha_KL   = alpha_y/(alpha_y+beta_y);   % = 0.3226 (preserve split K/L interne)
% m_Y calibre dans la boucle SS via FOC energie (Step C) ; sigma_Y->1 => m_Y -> 0.93
eta_X      = 0.60;       % export price elasticity (Henriet et al.)
theta      = 0.20;       % benchmark : 20% R_tau vers R&D
phi      = 0.80;       % benchmark : 80% transferts vers HtM

% Politique : taxe carbone initiale et finale
tau_ss     = 0.05;       % taxe carbone au SS initial (tau_init)
P_oil_ss   = 1.0;        % prix petrole normalise

% =========================================================================
% CALIBRATION ITERATIVE DU SS
% =========================================================================
fprintf('\n=== CALIBRATION V30-CES ===\n');

% --- Etape 1 : taux de rendement du capital ---
r_k_ss   = gz^sigma / beta - (1-delta);
fprintf('r_k_ss     = %.6f  (cible : r reelle ~1.5%%/trim)\n', r_k_ss);

% --- Etape 2 : normalisation output ---
Y_net_ss = 1.0;   % normalisation

% --- Etape 3 : facteur capital et investissement ---
K_ss     = alpha_y * Y_net_ss / r_k_ss;
I_ss     = (gz - (1-delta)) * K_ss;     % de l'accumulation : gz*K = (1-delta)*K + I
fprintf('K_ss       = %.6f\n', K_ss);
fprintf('I_ss       = %.6f  (I/Y = %.2f%%)\n', I_ss, 100*I_ss/Y_net_ss);

% --- Etape 4 : energie en production ---
% FOC energie : P_e*E_y = (1-alpha_y-beta_y)*Y_net
energy_share = 1 - alpha_y - beta_y;   % = 0.07
fprintf('energy_share = %.2f (= 1 - %.2f - %.2f) ✓ somme = 1\n', energy_share, alpha_y, beta_y);

% On a besoin du prix energie SS pour calibrer E_y_ss
% Prix energie propre : P_c = w_u/A_c, prix sale : P_d = P_oil + tau_ss
% Au SS sans friction mix (m constant), le CES donne P_e depuis P_c et P_d
% D'abord : taux de change er_ss (calibre D_ex_ss apres), poser er_ss = 1
er_ss    = 1.0;
P_d_ss   = er_ss * P_oil_ss + tau_ss;   % = 1.05

% On initialise une boucle iterative pour les marches du travail et energie
A_c_ss   = 1.0;   % normalisation technologie verte (cible)

% Guess initial
L_H_u_ss = 1/3;
L_R_s_ss = 1/3;

% Guess initial pour variables de convergence (AVANT la boucle)
R_tau_ss  = tau_ss * 0.004;
L_R_c_ss  = theta * R_tau_ss / 1.5;
L_u_c_ss  = 0.005;
T_H_ss    = phi*(1-theta)*R_tau_ss/lambda_htM;
T_R_ss    = (1-phi)*(1-theta)*R_tau_ss/(1-lambda_htM);

% Initialisation W_a_ss : approximation CD pour premier tour de boucle CES
W_a_ss = beta_y * Y_net_ss / 0.47;   % L_a ~ 0.47 (guess initial)
m_Y    = alpha_y + beta_y;            % = 0.93 (valeur CD initiale)

% --- Pre-calibration fonction dommages (cible fixe 1% GDP a +2°C, X=790 GtCO2) ---
% NOTE CRITIQUE (correction bug spec CES) : le document suppose Dam_ss=0 au SS
% initial, mais le modele reel a Dam_ss ~ 0.133 (X_ss~2881 GtCO2, stock de
% carbone atmospherique deja accumule). Il faut donc utiliser Y_g_ss =
% Y_net_ss/(1-Dam_ss) -- et non Y_net_ss -- dans les FOC de calibration CES
% (Steps C/D ci-dessous). Dam_ss depend de E_d_ss qui n'est connu qu'a la fin
% de CETTE iteration : on utilise donc Dam_ss_prev (retard d'une iteration),
% mis a jour a la fin de la boucle une fois E_d_ss calcule.
X_2deg   = 790;
d2       = 0.01 / X_2deg^2;
d1       = 0;
d0       = 0;
e_ROW_ss = 10;
Dam_ss_prev = 0.0;   % guess initial (converge vers Dam_ss au fil des iterations)
rho_Y = 1 - 1/sigma_Y;   % parametre CES (Y_g = ass*(m_Y*V^rho+(1-m_Y)*E^rho)^(1/rho))
ass_prev = 1.0;   % guess initial pour le facteur TFP ass^rho dans les FOC (converge en boucle)

for iter_main = 1:500
    % Travail total dans production
    L_s_y_ss  = (1-lambda_htM)*L_R_s_ss - L_R_c_ss;
    L_u_y_ss  = lambda_htM*L_H_u_ss - L_u_c_ss;
    if L_s_y_ss <= 0; L_s_y_ss = 0.01; end
    if L_u_y_ss <= 0; L_u_y_ss = 0.01; end

    % Agregat CES travail (CES sigma_L)
    % L_a = [omega_L*L_s_y^((sigma_L-1)/sigma_L) + (1-omega_L)*L_u_y^((sigma_L-1)/sigma_L)]^(sigma_L/(sigma_L-1))
    % Pour calibrer omega_L : cible prime salariale w_s/w_u = 1.60
    % FOC : w_s/w_u = (omega_L/(1-omega_L)) * (L_u_y/L_s_y)^(1/sigma_L)
    % => omega_L = (w_s/w_u) * (L_s_y/L_u_y)^(1/sigma_L) / (1 + (w_s/w_u)*(L_s_y/L_u_y)^(1/sigma_L))
    skill_premium_target = 1.60;
    ratio_Lu_Ls = L_u_y_ss / L_s_y_ss;
    omega_L    = skill_premium_target * ratio_Lu_Ls^(-1/sigma_L) / (1 + skill_premium_target * ratio_Lu_Ls^(-1/sigma_L));
    nu_L = (sigma_L-1)/sigma_L;
    L_a_ss = (omega_L*L_s_y_ss^nu_L + (1-omega_L)*L_u_y_ss^nu_L)^(1/nu_L);

    % Salaires (W_a_ss de l'iteration precedente -- mise a jour CES en aval)
    w_s_ss = W_a_ss * omega_L * (L_a_ss/L_s_y_ss)^(1/sigma_L);
    w_u_ss = W_a_ss * (1-omega_L) * (L_a_ss/L_u_y_ss)^(1/sigma_L);

    % Prix energie propre (V30 : Cobb-Douglas unit cost)
    P_c_ss = (1/A_c_ss) * (r_k_ss/alpha_c)^alpha_c * (w_u_ss/(1-alpha_c))^(1-alpha_c);
    % Prix CES energie (retailer SS sans friction)
    % Au SS : m constant => pas de terme friction
    % P_e FOC propre  : P_e*alpha_E*(E/E_c)^(1/rho) = P_c
    % P_e FOC sale    : P_e*(1-alpha_E)*(E/E_d)^(1/rho) = P_d
    % Formule prix minimum CES :
    P_e_ss  = (alpha_E^rho_E * P_c_ss^(1-rho_E) + (1-alpha_E)^rho_E * P_d_ss^(1-rho_E))^(1/(1-rho_E));

    % Energie en production (cible part 7%, preservee)
    E_y_ss  = energy_share * Y_net_ss / P_e_ss;

    % --- CES production finale Steps A-D ---
    % Step A : composite KL
    V_KL_ss = K_ss^alpha_KL * L_a_ss^(1-alpha_KL);
    % Y_g_ss approx (utilise Dam_ss retarde de l'iteration precedente ;
    % converge vers la valeur vraie car Dam_ss est mis a jour a chaque tour)
    Y_g_ss_approx = Y_net_ss / max(0.001, 1 - Dam_ss_prev);
    % Step C : m_Y via FOC energie (CES) -- MPK/MPE standard : dY_g/dX = a*ass^rho*(Y_g/X)^(1/sigma_Y)
    % => P_e = (1-m_Y)*ass_prev^rho_Y*(Y_g/E_y)^(1/sigma_Y)*(1-Dam)  (ass retarde, converge en boucle)
    m_Y = 1 - P_e_ss / ( ass_prev^rho_Y * (Y_g_ss_approx / E_y_ss)^(1/sigma_Y) * max(0.001, 1 - Dam_ss_prev) );
    m_Y = max(0.01, min(0.9999, m_Y));
    % Inner K_ss : point fixe FOC capital CES (r_k = alpha_KL*m_Y*ass^rho*(Y_g/V)^(1/sigma_Y)*(V/K)*(1-Dam))
    for iter_K = 1:100
        V_tmp = K_ss^alpha_KL * L_a_ss^(1-alpha_KL);
        K_new = alpha_KL * m_Y * ass_prev^rho_Y * (Y_g_ss_approx/V_tmp)^(1/sigma_Y) * V_tmp * (1-Dam_ss_prev) / r_k_ss;
        err_K = abs(K_new - K_ss);
        K_ss  = 0.5*K_ss + 0.5*K_new;
        if err_K < 1e-12; break; end
    end
    V_KL_ss = K_ss^alpha_KL * L_a_ss^(1-alpha_KL);
    I_ss    = (gz - (1-delta)) * K_ss;   % mise a jour avec K_y CES
    % Step D : salaire agrege CES (W_a = (1-alpha_KL)*m_Y*ass^rho*(Y_g/V)^(1/sigma_Y)*(V/L_a)*(1-Dam))
    W_a_ss_new  = (1-alpha_KL) * m_Y * ass_prev^rho_Y * (Y_g_ss_approx/V_KL_ss)^(1/sigma_Y) * V_KL_ss * (1-Dam_ss_prev) / L_a_ss;
    W_a_ss_prev = W_a_ss;
    W_a_ss      = 0.8*W_a_ss + 0.2*W_a_ss_new;   % mise a jour amortie

    % Energie propre et sale (allocation CES)
    % E_c/E_d = (alpha_E/(1-alpha_E))^rho_E * (P_d/P_c)^rho_E
    ratio_Ec_Ed = (alpha_E/(1-alpha_E))^rho_E * (P_d_ss/P_c_ss)^rho_E;
    % E_c = ratio_Ec_Ed * E_d
    % E = [alpha_E*E_c^nu_rho + (1-alpha_E)*E_d^nu_rho]^(1/nu_rho) where nu_rho = (rho-1)/rho
    % => E_d from E_y share (simplification SS : E ~= E_y + E_households)
    % Demande menages energie (Stone-Geary CES) -- bootstrap
    P_Z_ss  = ((1-omega) + omega*P_e_ss^(1-epsilon_E))^(1/(1-epsilon_E));

    % Transferts (utilise R_tau_ss de l'iteration precedente)
    T_H_ss     = phi * (1-theta) * R_tau_ss / lambda_htM;
    T_R_ss     = (1-phi) * (1-theta) * R_tau_ss / (1-lambda_htM);
    L_R_c_ss   = theta * R_tau_ss / w_s_ss;

    % Budget HtM -- cible L_H_u = 1/3 (chi calibre apres la boucle)
    ebar = 0.044;  % V30 : plancher subsistance energie (vs 0.05 en V29)
    L_H_u_ss = 1/3;
    budget_H  = w_u_ss*L_H_u_ss + T_H_ss - P_e_ss*ebar;
    if budget_H <= 0; ebar = 0.5*(w_u_ss*L_H_u_ss + T_H_ss)/P_e_ss; budget_H = w_u_ss*L_H_u_ss + T_H_ss - P_e_ss*ebar; end
    Z_H_ss = budget_H / P_Z_ss;

    % Budget Ricardien -- cible L_R_s = 1/3
    L_R_s_ss = 1/3;
    income_R  = w_s_ss*L_R_s_ss + r_k_ss*K_ss + T_R_ss - P_e_ss*ebar - I_ss;
    if income_R <= 0; income_R = 0.1; end
    Z_R_ss    = income_R / P_Z_ss;
    L_R_s_new = 1/3;

    % Energie menages (demand Stone-Geary)
    C_H_e_ss = ebar + omega*(P_Z_ss/P_e_ss)^epsilon_E * Z_H_ss;
    C_R_e_ss = ebar + omega*(P_Z_ss/P_e_ss)^epsilon_E * Z_R_ss;
    C_H_n_ss = (1-omega)*P_Z_ss^epsilon_E * Z_H_ss;
    C_R_n_ss = (1-omega)*P_Z_ss^epsilon_E * Z_R_ss;

    % Energie totale
    E_hh_ss = lambda_htM*C_H_e_ss + (1-lambda_htM)*C_R_e_ss;
    E_ss    = E_hh_ss + E_y_ss;
    % Decomposition E_c, E_d
    E_d_ss  = E_ss / (1 + ratio_Ec_Ed);
    E_c_ss  = ratio_Ec_Ed * E_d_ss;
    m_ss    = E_c_ss / E_d_ss;   % = ratio_Ec_Ed

    % Mise a jour Dam_ss_prev (pour le bloc CES de l'iteration suivante) :
    % bilan carbone coherent avec E_d_ss courant
    X_ss_iter   = (phi_d * E_d_ss + e_ROW_ss) / delta_clim;
    Dam_ss_prev = d0 + d1*X_ss_iter + d2*X_ss_iter^2;

    % Mise a jour ass_prev (facteur TFP^rho dans les FOC CES, point fixe amorti) :
    % ass = Y_g_ss / (m_Y*V_KL^rho + (1-m_Y)*E_y^rho)^(1/rho) (definition residuelle CES)
    Y_g_ss_iter = Y_net_ss / max(0.001, 1 - Dam_ss_prev);
    S_iter      = m_Y*V_KL_ss^rho_Y + (1-m_Y)*E_y_ss^rho_Y;
    ass_new     = Y_g_ss_iter / S_iter^(1/rho_Y);
    ass_prev    = 0.5*ass_prev + 0.5*ass_new;

    % Mise a jour L_u_c (V30 : demande facteur travail, pas fonction production directe)
    L_u_c_ss_new = (1-alpha_c) * P_c_ss * E_c_ss / w_u_ss;

    % Recettes taxe et allocations
    R_tau_ss_new = tau_ss * E_d_ss;
    T_H_ss_new   = phi*(1-theta)*R_tau_ss_new / lambda_htM;
    T_R_ss_new   = (1-phi)*(1-theta)*R_tau_ss_new / (1-lambda_htM);
    L_R_c_ss_new = theta * R_tau_ss_new / w_s_ss;

    % Convergence
    err_main = max(abs([L_u_c_ss_new - L_u_c_ss, R_tau_ss_new - R_tau_ss, ...
                        L_R_c_ss_new - L_R_c_ss, W_a_ss_new - W_a_ss_prev]));

    % L_R_s_ss reste fixee a 1/3 (cible)
    L_u_c_ss = 0.5*L_u_c_ss + 0.5*L_u_c_ss_new;
    R_tau_ss = 0.5*R_tau_ss + 0.5*R_tau_ss_new;
    L_R_c_ss = 0.5*L_R_c_ss + 0.5*L_R_c_ss_new;
    T_H_ss   = T_H_ss_new;
    T_R_ss   = T_R_ss_new;

    if err_main < 1e-9; fprintf('Convergence SS atteinte en %d iterations\n', iter_main); break; end
end

% --- V30 : Capital dans energie propre et total ---
K_c_ss     = alpha_c * P_c_ss * E_c_ss / r_k_ss;
K_total_ss = K_ss + K_c_ss;            % K_ss = K_y (bien final seulement)
I_ss       = (gz - (1-delta)) * K_total_ss;   % investissement couvre K_y + K_c
% Budget Ricardien corriger (revenu capital = r_k * K_total)
income_R_corr = w_s_ss*L_R_s_ss + r_k_ss*K_total_ss + T_R_ss - P_e_ss*ebar - I_ss;
if income_R_corr > 0
    Z_R_ss = income_R_corr / P_Z_ss;
end
fprintf('K_c_ss     = %.6f (V30: capital energie propre)\n', K_c_ss);
fprintf('K_total_ss = %.6f (K_y + K_c)\n', K_total_ss);
fprintf('I_ss (V30) = %.6f (couvre K_y + K_c)\n', I_ss);

% --- Calibration de chi (poids travail, cible L_ss = 1/3 pour chaque type) ---
% HtM labor FOC : chi*L_H_u^psi = Z_H^(-sigma)*w_u/P_Z
chi_H = Z_H_ss^(-sigma)*w_u_ss / (P_Z_ss * L_H_u_ss^psi);
% Ricardian labor FOC : chi*L_R_s^psi = Z_R^(-sigma)*w_s/P_Z
chi_R = Z_R_ss^(-sigma)*w_s_ss / (P_Z_ss * L_R_s_ss^psi);
chi   = 0.5*(chi_H + chi_R);  % moyenne (devrait etre proche si bonne calibration)
fprintf('chi (H)    = %.6f | chi (R) = %.6f | chi retenu = %.6f\n', chi_H, chi_R, chi);

% --- Calibration de ebar (cible burden HtM ~19%) ---
% burden_H = P_e*C_H_e / (C_H_n + P_e*C_H_e)
% On ajuste ebar jusqu'a atteindre burden_H ~ 0.19
burden_H_actual = P_e_ss*C_H_e_ss / (C_H_n_ss + P_e_ss*C_H_e_ss);
burden_R_actual = P_e_ss*C_R_e_ss / (C_R_n_ss + P_e_ss*C_R_e_ss);
regressivity    = burden_H_actual / burden_R_actual;
fprintf('burden_H   = %.4f (cible ~0.19)\n', burden_H_actual);
fprintf('burden_R   = %.4f (cible ~0.07)\n', burden_R_actual);
fprintf('regressivite = %.4f (cible ~2.5)\n', regressivity);
fprintf('[Note: ajuster ebar manuellement si burden_H trop loin de 0.19]\n');

% --- Calibration eta_c (SS innovation verte : delta_c*A_c_ss = eta_c*L_R_c_ss) ---
eta_c = delta_c * A_c_ss / L_R_c_ss;
fprintf('eta_c      = %.6f (calibre EE : delta_c*A_c_ss/L_R_c_ss)\n', eta_c);

% --- Calibration damage function (DOIT etre fait AVANT ass) ---
% Cible : 1% GDP a +2 deg C (X = 790 GtCO2)
X_2deg = 790;
d2 = 0.01 / X_2deg^2;   % approx 1.603e-8
d1 = 0;
d0 = 0;
fprintf('Dam a +2°C = %.4f (cible 1%% GDP, X=790 GtCO2)\n', d0+d1*X_2deg+d2*X_2deg^2);

% --- Calibration ass (TFP neutre) ---
% X_ss consistant avec le bilan carbone du modele : delta_clim*X = e_em + e_ROW
e_ROW_ss = 10;
X_ss = (phi_d * E_d_ss + e_ROW_ss) / delta_clim;  % ~2891 GtCO2
Dam_ss = d0 + d1*X_ss + d2*X_ss^2;
Y_g_ss = Y_net_ss / (1 - Dam_ss);
% Step E : calibrer ass (TFP neutre, formule CES)
V_KL_ss = K_ss^alpha_KL * L_a_ss^(1-alpha_KL);   % valeur finale apres convergence
ass = Y_g_ss / ( m_Y*V_KL_ss^(1-1/sigma_Y) + (1-m_Y)*E_y_ss^(1-1/sigma_Y) )^(1/(1-1/sigma_Y));
fprintf('ass        = %.6f (TFP CES normalisee)\n', ass);
% Step F : verification ex post FOC capital (ecart < 1e-6 requis)
r_k_check = alpha_KL * m_Y * ass^(1-1/sigma_Y) * (Y_g_ss/V_KL_ss)^(1/sigma_Y) * V_KL_ss*(1-Dam_ss)/K_ss;
fprintf('r_k FOC check : r_k_ss=%.6f vs CES=%.6f (ecart=%.2e)\n', r_k_ss, r_k_check, abs(r_k_ss-r_k_check));
fprintf('m_Y        = %.6f | alpha_KL = %.6f | sigma_Y = %.4f\n', m_Y, alpha_KL, sigma_Y);
fprintf('V_KL_ss    = %.6f | K_ss CES = %.6f (vs K_ss CD = %.6f)\n', V_KL_ss, K_ss, alpha_y*Y_net_ss/r_k_ss);
fprintf('Dam_ss     = %.6f (Dam SS, X_ss=%.1f GtCO2)\n', Dam_ss, X_ss);
fprintf('X_ss       = %.1f (SS consistant avec bilan carbone)\n', X_ss);

% --- Calibration D_ex_ss (demande etrangere, cible er_ss=1) ---
D_ex_ss = P_oil_ss * er_ss * E_d_ss * er_ss^eta_X;
fprintf('D_ex_ss    = %.6f (exportations au SS = facture petrole)\n', D_ex_ss);

% --- Welfare SS (pour CEV) ---
% U_H_ss = Z_H_ss^(1-sigma)/(1-sigma) - chi*L_H_u_ss^(1+psi)/(1+psi)
U_H_ss = Z_H_ss^(1-sigma)/(1-sigma) - chi*L_H_u_ss^(1+psi)/(1+psi);
U_R_ss = Z_R_ss^(1-sigma)/(1-sigma) - chi*L_R_s_ss^(1+psi)/(1+psi);
% W_ss = U_ss / (1 - beta*gz^(1-sigma))
disc_W = beta * gz^(1-sigma);
if abs(1-disc_W) < 1e-8; error('Facteur actualisation welfare = 1, instable'); end
W_H_ss = U_H_ss / (1 - disc_W);
W_R_ss = U_R_ss / (1 - disc_W);
fprintf('W_H_ss     = %.4f | W_R_ss = %.4f\n', W_H_ss, W_R_ss);
fprintf('Facteur actualisation welfare : beta*gz^(1-sigma) = %.6f\n', disc_W);

% --- Taux reel SS (pour initval) ---
rr_ss = gz^sigma / beta;  % taux reel brut de court terme
r_ss  = rr_ss;            % modele reel (pas d'inflation)

% --- Agregats SS ---
C_n_ss = lambda_htM*C_H_n_ss + (1-lambda_htM)*C_R_n_ss;
EX_ss  = P_oil_ss * er_ss * E_d_ss;
Lambda_R_ss = Z_R_ss^(-sigma) / P_Z_ss;

% --- Tobin q SS ---
q_K_ss = 1 / (1 - phi_k*(I_ss/K_total_ss - delta));   % V30 : K_total
gamma_agg_ss = gz;

fprintf('\n--- RESUME SS V30-CES ---\n');
fprintf('Y_net_ss   = %.4f | K_ss     = %.4f | I_ss    = %.4f\n', Y_net_ss, K_ss, I_ss);
fprintf('w_u_ss     = %.4f | w_s_ss   = %.4f | prem.sal = %.4f\n', w_u_ss, w_s_ss, w_s_ss/w_u_ss);
fprintf('L_H_u_ss   = %.4f | L_R_s_ss = %.4f\n', L_H_u_ss, L_R_s_ss);
fprintf('P_e_ss     = %.4f | P_c_ss   = %.4f | P_d_ss  = %.4f\n', P_e_ss, P_c_ss, P_d_ss);
fprintf('E_ss       = %.6f | E_c_ss   = %.6f | E_d_ss  = %.6f\n', E_ss, E_c_ss, E_d_ss);
fprintf('A_c_ss     = %.4f | L_R_c_ss = %.6f | L_u_c_ss= %.6f\n', A_c_ss, L_R_c_ss, L_u_c_ss);
fprintf('R_tau_ss   = %.6f | T_H_ss   = %.6f | T_R_ss  = %.6f\n', R_tau_ss, T_H_ss, T_R_ss);
fprintf('er_ss      = %.4f | EX_ss    = %.6f\n', er_ss, EX_ss);
fprintf('P_Z_ss     = %.4f | Z_H_ss   = %.4f | Z_R_ss  = %.4f\n', P_Z_ss, Z_H_ss, Z_R_ss);
fprintf('omega_L    = %.4f (calibre prime salariale)\n', omega_L);
fprintf('chi        = %.4f (calibre L_ss~1/3)\n', chi);
fprintf('eta_c      = %.6f (calibre A_c_ss=1)\n', eta_c);
fprintf('D_ex_ss    = %.6f (calibre er_ss=1)\n', D_ex_ss);
fprintf('=====================\n');

%%
%%%%%%%%%%%%% MODEL BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model;
% =========================================================================
% BLOC 1 : PREFERENCES STONE-GEARY CES
% =========================================================================

% [1] Indice de prix CES composite energie/non-energie
% P_Z = [(1-omega) + omega*P_e^(1-epsilon_E)]^(1/(1-epsilon_E))
P_Z = ((1-omega) + omega*P_e^(1-epsilon_E))^(1/(1-epsilon_E));

% [2] Demande non-energie HtM (numeraire P_n = 1)
C_H_n = (1-omega)*P_Z^epsilon_E * Z_H;

% [3] Demande energie HtM (supranumeraire + plancher)
% C_H_e = ebar + omega*(P_Z/P_e)^epsilon_E * Z_H
C_H_e = ebar + omega*(P_Z/P_e)^epsilon_E * Z_H;

% [4] Demande non-energie Ricardien
C_R_n = (1-omega)*P_Z^epsilon_E * Z_R;

% [5] Demande energie Ricardien
C_R_e = ebar + omega*(P_Z/P_e)^epsilon_E * Z_R;

% [6] Budget HtM (en termes Z_H)
% P_Z*Z_H + P_e*ebar = w_u*L_H_u + T_H
% => P_Z*Z_H = w_u*L_H_u + T_H - P_e*ebar
P_Z*Z_H = w_u*L_H_u + T_H - P_e*ebar;

% [7] FOC travail HtM (Stone-Geary : MU conso = Z_H^{-sigma}/P_Z)
% chi*L_H_u^psi = Z_H^(-sigma)*w_u/P_Z
chi*L_H_u^psi = Z_H^(-sigma)*w_u/P_Z;

% [8] Utilite marginale revenu Ricardien
Lambda_R = Z_R^(-sigma)/P_Z;

% [9] Budget Ricardien -- V30 : capital income = r_k * K_total = r_k * K(-1)
% (K(-1) est le capital total fin de periode t-1 = K_y + K_c, r_k commun aux deux secteurs)
P_Z*Z_R = w_s*L_R_s + r_k*K(-1) + T_R - P_e*ebar - I;

% [10] FOC travail Ricardien
chi*L_R_s^psi = Z_R^(-sigma)*w_s/P_Z;

% [11] q de Tobin (condition d'optimalite investissement, couts I/K)
% q_K*(1 - phi_k*(I/K(-1) - delta)) = 1
q_K*(1 - phi_k*(I/K(-1) - delta)) = 1;

% [12] Accumulation du capital (detrende par gz)
% gz*K = (1-delta)*K(-1) + I - phi_k/2*(I/K(-1) - delta)^2 * K(-1)
gz*K = (1-delta)*K(-1) + I - phi_k/2*(I/K(-1) - delta)^2*K(-1);

% [13] Equation d'Euler (investissement + Tobin Q forward-looking)
% Lambda_R*q_K = beta*gz^(-sigma)*Lambda_R(+1)*(r_k(+1) + q_K(+1)*(1-delta))
Lambda_R*q_K = beta*gz^(-sigma)*Lambda_R(+1)*(r_k(+1) + q_K(+1)*(1-delta));

% =========================================================================
% BLOC 2 : RETAILER ENERGIE (CES + friction de mix)
% =========================================================================

% [14] Aggregateur CES energie
% E = [alpha_E*E_c^((rho_E-1)/rho_E) + (1-alpha_E)*E_d^((rho_E-1)/rho_E)]^(rho_E/(rho_E-1))
E = (alpha_E*E_c^((rho_E-1)/rho_E) + (1-alpha_E)*E_d^((rho_E-1)/rho_E))^(rho_E/(rho_E-1));

% [15] Ratio de mix energie (variable d'etat : m(-1) dans les FOC)
m = E_c/E_d;

% [16] Prix energie sale (taxee)
P_d = er*P_oil + tau;

% [17] Prix energie propre -- V30 : Cobb-Douglas unit cost (vs w_u/A_c en V29)
% P_c = (1/A_c)*(r_k/alpha_c)^alpha_c * (w_u/(1-alpha_c))^(1-alpha_c)
P_c = (1/A_c) * (r_k/alpha_c)^alpha_c * (w_u/(1-alpha_c))^(1-alpha_c);

% [18] FOC retailer pour energie propre (avec friction mix)
% P_e*alpha_E*(E/E_c)^(1/rho_E) = P_c + phi_mix*(m - m(-1))*P_e*E/E_d
P_e*alpha_E*(E/E_c)^(1/rho_E) = P_c + phi_mix*(m - m(-1))*P_e*E/E_d;

% [19] FOC retailer pour energie sale (avec friction mix)
% P_e*(1-alpha_E)*(E/E_d)^(1/rho_E) = P_d - phi_mix*(m - m(-1))*E_c*P_e*E/E_d^2
P_e*(1-alpha_E)*(E/E_d)^(1/rho_E) = P_d - phi_mix*(m - m(-1))*E_c*P_e*E/E_d^2;

% =========================================================================
% BLOC 3 : PRODUCTION FINALE (KLE, parts somment a 1)
% =========================================================================

% [20] Composite capital-travail (Cobb-Douglas interne, rendements constants)
V_KL = (K(-1)-K_c)^alpha_KL * L_a^(1-alpha_KL);

% [20b] Production finale CES gigogne : composite KL + energie
Y_g = ass * ( m_Y*V_KL^(1-1/sigma_Y) + (1-m_Y)*E_y^(1-1/sigma_Y) )^(1/(1-1/sigma_Y));

% [21] Output net de dommages (inchange)
Y_net = (1 - Dam) * Y_g;

% [22] FOC capital (CES externe, via composite V_KL)
% NOTE : facteur ass^(1-1/sigma_Y) requis par la derivee CES standard
% MPK = m_Y*ass^rho*(Y_g/V)^(1/sigma_Y) ; le second terme est V_KL/K (pas Y_g/K)
r_k = alpha_KL * m_Y * ass^(1-1/sigma_Y) * (Y_g/V_KL)^(1/sigma_Y) * V_KL/(K(-1)-K_c) * (1-Dam);

% [23] FOC energie en production (CES externe)
P_e = (1-m_Y) * ass^(1-1/sigma_Y) * (Y_g/E_y)^(1/sigma_Y) * (1-Dam);

% [24] FOC travail agrege (CES externe, via composite V_KL)
W_a = (1-alpha_KL) * m_Y * ass^(1-1/sigma_Y) * (Y_g/V_KL)^(1/sigma_Y) * V_KL/L_a * (1-Dam);

% [25] Agregat CES de travail (sigma_L)
L_a = (omega_L*L_s_y^((sigma_L-1)/sigma_L) + (1-omega_L)*L_u_y^((sigma_L-1)/sigma_L))^(sigma_L/(sigma_L-1));

% [26] FOC travail qualifie en production
w_s = W_a*omega_L*(L_a/L_s_y)^(1/sigma_L);

% [27] FOC travail non-qualifie en production
w_u = W_a*(1-omega_L)*(L_a/L_u_y)^(1/sigma_L);

% [28] Demande facteur travail non-qualifie dans energie propre -- V30
% L_u_c = (1-alpha_c)*P_c*E_c/w_u  (part travail du cout, profit nul CRS)
L_u_c = (1-alpha_c)*P_c*E_c/w_u;

% [46] Demande facteur capital dans energie propre -- V30 (nouvelle equation)
% K_c = alpha_c*P_c*E_c/r_k  (part capital du cout, profit nul CRS)
K_c = alpha_c*P_c*E_c/r_k;

% =========================================================================
% BLOC 4 : INNOVATION VERTE ET CROISSANCE
% =========================================================================

% [29] Stock de technologie verte (depreciating knowledge stock)
% A_c = (1-delta_c)*A_c(-1) + eta_c*L_R_c
A_c = (1-delta_c)*A_c(-1) + eta_c*L_R_c;

% [30] Taux de croissance semi-endogene (reporting seulement)
gamma_agg = gz*Y_net/Y_net(-1);

% =========================================================================
% BLOC 5 : CLIMAT (DICE quadratique)
% =========================================================================

% [31] Emissions domestiques
e_em = phi_d*E_d;

% [32] Stock de carbone atmospherique
X = (1-delta_clim)*X(-1) + e_em + e_ROW;

% [33] Fonction de dommages DICE
Dam = d0 + d1*X + d2*X^2;

% =========================================================================
% BLOC 6 : GOUVERNEMENT -- Taxe carbone et recyclage
% =========================================================================
% Budget : theta*R_tau -> R&D vert
%          phi*(1-theta)*R_tau -> transferts HtM
%          (1-phi)*(1-theta)*R_tau -> transferts Ricardiens
% Verification : theta + phi*(1-theta) + (1-phi)*(1-theta) = 1 ✓

% [34] Recettes taxe carbone
R_tau = tau*E_d;

% [35] Subvention R&D vert (financement L_R_c)
L_R_c = theta*R_tau/w_s;

% [36] Transferts vers menages HtM (per capita)
T_H = phi*(1-theta)*R_tau/lambda_htM;

% [37] Transferts vers Ricardiens (per capita)
T_R = (1-phi)*(1-theta)*R_tau/(1-lambda_htM);

% =========================================================================
% BLOC 7 : EQUILIBRES DES MARCHES + OPEN ECONOMY
% =========================================================================

% [38] Agregat conso non-energie
C_n = lambda_htM*C_H_n + (1-lambda_htM)*C_R_n;

% Note : Marche des biens (Y_net = C_n + I + adj_cost) implicite par Walras.
% Verification post-simulation : residu doit etre < 1e-6.

% [39] Marche du travail non-qualifie
% lambda*L_H_u = L_u_y + L_u_c
lambda_htM*L_H_u = L_u_y + L_u_c;

% [40] Marche du travail qualifie
% (1-lambda)*L_R_s = L_s_y + L_R_c
(1-lambda_htM)*L_R_s = L_s_y + L_R_c;

% [41] Marche de l'energie (demande totale = conso menages + production)
% E = lambda*C_H_e + (1-lambda)*C_R_e + E_y
E = lambda_htM*C_H_e + (1-lambda_htM)*C_R_e + E_y;

% [42] Demande etrangere d'exportations
EX = D_ex_ss*er^(-eta_X);

% [43] Fermeture open economy (exports = facture petrole, trade balance = 0)
% EX = P_oil*er*E_d => financement importe par exports biens finaux
EX = P_oil*er*E_d;

% =========================================================================
% BLOC 8 : BIEN-ETRE (recursions, pour CEV)
% =========================================================================

% [44] Bien-etre HtM (recursion backward de la valeur actuelle)
% W_H = [Z_H^(1-sigma)/(1-sigma) - chi*L_H_u^(1+psi)/(1+psi)] + beta*gz^(1-sigma)*W_H(+1)
W_H = Z_H^(1-sigma)/(1-sigma) - chi*L_H_u^(1+psi)/(1+psi) + beta*gz^(1-sigma)*W_H(+1);

% [45] Bien-etre Ricardien
W_R = Z_R^(1-sigma)/(1-sigma) - chi*L_R_s^(1+psi)/(1+psi) + beta*gz^(1-sigma)*W_R(+1);

end;

%%
%%%%%%%%%%%%% SIMULATION PERFECT FORESIGHT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Horizon : T_sim = 200 trimestres (50 ans)
% Taxe carbone : rampe convexe de tau_init=0.05 a tau_final=0.40 sur 120 trim. (30 ans)

@#define T_sim = 200

% Conditions initiales (SS initial, tau=0.05)
initval;
tau      = tau_ss;
P_oil    = P_oil_ss;
e_ROW    = e_ROW_ss;
er       = er_ss;
P_d      = er_ss*P_oil_ss + tau_ss;
A_c      = A_c_ss;
P_c      = P_c_ss;
P_e      = P_e_ss;
P_Z      = P_Z_ss;
m        = m_ss;
K        = K_total_ss;   % V30 : K_total = K_y + K_c
K_c      = K_c_ss;       % V30 : capital energie propre SS initial
I        = I_ss;
r_k      = r_k_ss;
q_K      = q_K_ss;
w_u      = w_u_ss;
w_s      = w_s_ss;
W_a      = W_a_ss;
L_a      = L_a_ss;
L_s_y    = L_s_y_ss;
L_u_y    = L_u_y_ss;
L_u_c    = L_u_c_ss;
L_H_u    = L_H_u_ss;
L_R_s    = L_R_s_ss;
Z_H      = Z_H_ss;
Z_R      = Z_R_ss;
C_H_n    = C_H_n_ss;
C_R_n    = C_R_n_ss;
C_H_e    = C_H_e_ss;
C_R_e    = C_R_e_ss;
Lambda_R = Lambda_R_ss;
C_n      = C_n_ss;
E        = E_ss;
E_c      = E_c_ss;
E_d      = E_d_ss;
E_y      = E_y_ss;
Y_g      = Y_g_ss;
Y_net    = Y_net_ss;
V_KL     = V_KL_ss;
L_R_c    = L_R_c_ss;
gamma_agg = gz_ss;
e_em     = phi_d*E_d_ss;
X        = X_ss;
Dam      = Dam_ss;
R_tau    = R_tau_ss;
T_H      = T_H_ss;
T_R      = T_R_ss;
EX       = EX_ss;
W_H      = W_H_ss;
W_R      = W_R_ss;
end;

% Conditions terminales (SS final, tau=0.40 par defaut -- surchargeable par
% les scripts de grille via regexprep sur tau_final_v, cf. plus bas)
endval;
tau      = 0.40;
P_oil    = P_oil_ss;
e_ROW    = e_ROW_ss;
end;

% Choc : passage de tau_init a tau_final sur 120 trimestres (rampe convexe)
% La trajectoire exacte est injectee plus bas, directement dans ce .mod,
% dans oo_.exo_simul apres perfect_foresight_setup (cf. bloc kappa_v=2)
shocks;
var tau;
periods 1;
values 0.05;
end;
% Note : le bloc ci-dessous (apres perfect_foresight_setup) surcharge
% oo_.exo_simul avec la trajectoire exacte (rampe convexe kappa=2)

%%
%%%%%%%%%%%%% SS NUMERIQUE (apres initval/endval/shocks) %%%%%%%%%%%%%%%%%%
% steady doit etre appele APRES initval et endval pour utiliser ces valeurs
% comme point de depart de Newton
steady;
check;

%%
%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perfect_foresight_setup(periods = @{T_sim});

%% Injecter la trajectoire LISSE de tau (rampe quadratique) dans oo_.exo_simul
% Remplace le saut abrupt provenant du bloc shocks (qui ne fixait que t=1)
tau_idx_v   = find(strcmp(M_.exo_names, 'tau'));
tau_init_v  = 0.05;
tau_final_v = 0.40;
T_carbon_v  = 120;
kappa_v     = 2;
T_sim_v     = @{T_sim};   % = 200
% oo_.exo_simul est (T_sim+2) x n_exo :
%   ligne 1 = t=0 (initval), lignes 2..T+1 = transition, ligne T+2 = endval
oo_.exo_simul(1, tau_idx_v) = tau_init_v;
for t_v = 1:T_sim_v
    if t_v <= T_carbon_v
        oo_.exo_simul(t_v+1, tau_idx_v) = tau_init_v + ...
            (tau_final_v - tau_init_v) * (t_v / T_carbon_v)^kappa_v;
    else
        oo_.exo_simul(t_v+1, tau_idx_v) = tau_final_v;
    end
end
oo_.exo_simul(T_sim_v+2, tau_idx_v) = tau_final_v;

perfect_foresight_solver(maxit = 100, tolf = 1e-8);

%%
%%%%%%%%%%%%% EXTRACTION RESULTATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if oo_.deterministic_simulation.status == 1

    T_sim_read = @{T_sim};
    get_idx = @(name) find(strcmp(M_.endo_names, name));

    t_q = (1:T_sim_read)';

    % --- Macro ---
    Y_net_sim  = oo_.endo_simul(get_idx('Y_net'),  2:T_sim_read+1)';
    C_n_sim    = oo_.endo_simul(get_idx('C_n'),    2:T_sim_read+1)';
    I_sim      = oo_.endo_simul(get_idx('I'),      2:T_sim_read+1)';
    w_u_sim    = oo_.endo_simul(get_idx('w_u'),    2:T_sim_read+1)';
    w_s_sim    = oo_.endo_simul(get_idx('w_s'),    2:T_sim_read+1)';

    % --- Energie ---
    P_e_sim    = oo_.endo_simul(get_idx('P_e'),    2:T_sim_read+1)';
    E_d_sim    = oo_.endo_simul(get_idx('E_d'),    2:T_sim_read+1)';
    E_c_sim    = oo_.endo_simul(get_idx('E_c'),    2:T_sim_read+1)';
    E_sim      = oo_.endo_simul(get_idx('E'),      2:T_sim_read+1)';
    m_sim      = oo_.endo_simul(get_idx('m'),      2:T_sim_read+1)';

    % --- Distribution et bien-etre ---
    C_H_e_sim  = oo_.endo_simul(get_idx('C_H_e'), 2:T_sim_read+1)';
    C_R_e_sim  = oo_.endo_simul(get_idx('C_R_e'), 2:T_sim_read+1)';
    C_H_n_sim  = oo_.endo_simul(get_idx('C_H_n'), 2:T_sim_read+1)';
    C_R_n_sim  = oo_.endo_simul(get_idx('C_R_n'), 2:T_sim_read+1)';
    Z_H_sim    = oo_.endo_simul(get_idx('Z_H'),   2:T_sim_read+1)';
    Z_R_sim    = oo_.endo_simul(get_idx('Z_R'),   2:T_sim_read+1)';
    W_H_sim    = oo_.endo_simul(get_idx('W_H'),   2:T_sim_read+1)';
    W_R_sim    = oo_.endo_simul(get_idx('W_R'),   2:T_sim_read+1)';
    T_H_sim    = oo_.endo_simul(get_idx('T_H'),   2:T_sim_read+1)';
    T_R_sim    = oo_.endo_simul(get_idx('T_R'),   2:T_sim_read+1)';

    % --- Innovation et emplois verts ---
    A_c_sim    = oo_.endo_simul(get_idx('A_c'),      2:T_sim_read+1)';
    L_R_c_sim  = oo_.endo_simul(get_idx('L_R_c'),    2:T_sim_read+1)';
    L_u_c_sim  = oo_.endo_simul(get_idx('L_u_c'),    2:T_sim_read+1)';
    gamma_sim  = oo_.endo_simul(get_idx('gamma_agg'),2:T_sim_read+1)';

    % --- Climat ---
    X_sim      = oo_.endo_simul(get_idx('X'),      2:T_sim_read+1)';
    Dam_sim    = oo_.endo_simul(get_idx('Dam'),    2:T_sim_read+1)';
    e_em_sim   = oo_.endo_simul(get_idx('e_em'),   2:T_sim_read+1)';

    % --- Burdens et regressivite ---
    burden_H_sim = P_e_sim.*C_H_e_sim ./ (C_H_n_sim + P_e_sim.*C_H_e_sim);
    burden_R_sim = P_e_sim.*C_R_e_sim ./ (C_R_n_sim + P_e_sim.*C_R_e_sim);
    regress_sim  = burden_H_sim ./ burden_R_sim;
    prem_sal_sim = w_s_sim ./ w_u_sim;

    % SS (pour normalisation)
    ss = oo_.steady_state;
    Y_ss   = ss(get_idx('Y_net'));
    Pe_ss  = ss(get_idx('P_e'));
    Ac_ss  = ss(get_idx('A_c'));
    WH_ss  = ss(get_idx('W_H'));
    WR_ss  = ss(get_idx('W_R'));

    fprintf('\n=== V30-CES RESULTATS A 30 ANS (periode 120) ===\n');
    fprintf('Y_net     : %+.2f%% vs SS\n', 100*(Y_net_sim(120)/Y_ss-1));
    fprintf('P_e       : %+.2f%% vs SS\n', 100*(P_e_sim(120)/Pe_ss-1));
    fprintf('A_c       : %+.2f%% vs SS\n', 100*(A_c_sim(120)/Ac_ss-1));
    fprintf('burden_H  : %.4f (SS: %.4f)\n', burden_H_sim(120), P_e_ss*C_H_e_ss/(C_H_n_ss+P_e_ss*C_H_e_ss));
    fprintf('regressiv : %.4f (SS: %.4f)\n', regress_sim(120), regressivity);

    % S-CEV (Supranumerary Consumption Equivalent Variation)
    % W^transition = W^ss * (1+mu)^(1-sigma) / (1-sigma) * (1-sigma)
    % Pour sigma != 1 : mu = (W_t1/W_ss)^(1/(1-sigma)) - 1
    % Attention : avec sigma=1.5 > 1, W < 0 => adapter le signe
    WH_t1 = W_H_sim(1);   WR_t1 = W_R_sim(1);
    if WH_ss < 0 && WH_t1 < 0
        % welfare negatif : ratio |W_t1/W_ss| eleve = mieux
        SCEV_H = 100*(abs(WH_t1/WH_ss)^(1/(1-sigma)) - 1);
        if WH_t1 > WH_ss; SCEV_H = abs(SCEV_H); else; SCEV_H = -abs(SCEV_H); end
    else
        SCEV_H = 100*((WH_t1/WH_ss)^(1/(1-sigma)) - 1);
    end
    if WR_ss < 0 && WR_t1 < 0
        SCEV_R = 100*(abs(WR_t1/WR_ss)^(1/(1-sigma)) - 1);
        if WR_t1 > WR_ss; SCEV_R = abs(SCEV_R); else; SCEV_R = -abs(SCEV_R); end
    else
        SCEV_R = 100*((WR_t1/WR_ss)^(1/(1-sigma)) - 1);
    end
    fprintf('S-CEV HtM : %+.4f%% | S-CEV Ric : %+.4f%%\n', SCEV_H, SCEV_R);

else
    warning('Simulation V30-CES n''a pas converge. Verifier initval/endval et les residus.');
end
