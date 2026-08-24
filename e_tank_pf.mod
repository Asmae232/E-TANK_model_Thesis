%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% E-TANK PERFECT FORESIGHT -- V14 ISSUE 2 (full steady-state consistency)
%%%% Author: Asmae EL MOUHSSINE -- June 2026
%%%%
%%%% FINAL PRODUCTION STRUCTURE (nested CES, Schubert's suggestion) :
%%%%   V_KL  = (K_y)^alpha_KL * L_a^(1-alpha_KL)                    (K-L composite, inner Cobb-Douglas)
%%%%   Y_g   = ass * [ m_Y*V_KL^(1-1/sigma_Y) + (1-m_Y)*E_y^(1-1/sigma_Y) ]^(1/(1-1/sigma_Y))
%%%%   Cobb-Douglas (historical benchmark) is the limiting case sigma_Y -> 1,
%%%%   with m_Y = alpha_y+beta_y = 0.93 and alpha_KL = alpha_y/(alpha_y+beta_y) = 0.3226
%%%%   (validated by Test 1: deviations < 0.01pp vs CD at sigma_Y=0.9999).
%%%%
%%%% DIFFERENCES V30 vs V29 (4 surgical changes, inherited from the CD base) :
%%%%   (1) Clean energy : E_c = A_c * K_c^alpha_c * L_u_c^(1-alpha_c)
%%%%       alpha_c = 0.20 -- corrects the skill-premium inversion
%%%%       P_c = Cobb-Douglas unit cost ; K = K_y + K_c (capital equilibrium)
%%%%   (2) sigma_L = 1.40 (vs 0.80) -- Katz-Murphy / Acemoglu-Autor
%%%%   (3) omega   = 0.012 -- fixed/exogenous preference parameter
%%%%   (4) ebar    = 0.044 (vs 0.05) -- subsistence floor
%%%%
%%%% ENDOGENOUS VARIABLES : 47 (46 base V30 + V_KL, nested CES composite)
%%%% EQUATIONS : 47 (goods market implicit by Walras)
%%%% VAREXO : tau, P_oil, e_ROW, theta_policy, phi_policy
%%%%
%%%% GENERAL EQUILIBRIUM :
%%%%   - Two household types : HtM (lambda=0.20) + Ricardians (1-lambda=0.80)
%%%%   - Nested CES final production (K-L)/Energy; SS shares calibrated to match CD
%%%%   - Carbon-tax recycling : R_tau = tau*E_d
%%%%       theta*R_tau -> green R&D subsidy
%%%%       phi*(1-theta)*R_tau -> HtM transfers
%%%%       (1-phi)*(1-theta)*R_tau -> Ricardian transfers
%%%%   - Open-economy closure : exports = oil bill (tb=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; warning off;

@#ifndef THETA_POLICY
  @#define THETA_POLICY = 0.20
@#endif
@#ifndef PHI_POLICY
  @#define PHI_POLICY = 0.80
@#endif
@#ifndef P_OIL_FINAL
  @#define P_OIL_FINAL = 1.0
@#endif
@#ifndef IDENT_INSTANT_SHOCK
  @#define IDENT_INSTANT_SHOCK = 0
@#endif
@#ifndef P_OIL_PATH_MODE
  @#define P_OIL_PATH_MODE = 0
@#endif
@#ifndef P_OIL_PERSISTENCE
  @#define P_OIL_PERSISTENCE = 0.95
@#endif

%%
%%%%%%%%%%%%% Endogenous variables (45) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var
% Preferences Stone-Geary CES (3)
P_Z         % CES price index energy/non-energy
Z_H         % HtM welfare composite (Stone-Geary CES)
Z_R         % Ricardian welfare composite

% Consumption by type and good (4)
C_H_n       % HtM non-energy consumption
C_R_n       % Ricardian non-energy consumption
C_H_e       % HtM energy consumption (= ebar + supranumeraire)
C_R_e       % Ricardian energy consumption

% Labor supply (2)
L_H_u       % HtM hours worked (unskilled)
L_R_s       % Ricardian hours worked (skilled)

% Ricardian decisions: capital, investment, MU (4)
Lambda_R    % Ricardian marginal utility of income
I           % investment (detrended)
K           % capital (detrended, end-of-period stock)
q_K         % Tobin's q (Ricardian)

% Energy prices and quantities (retailer CES + mix friction) (9)
E           % total energy (CES composite)
E_c         % clean energy (domestic)
E_d         % dirty energy (imported)
E_y         % energy used in production
m           % mix ratio E_c/E_d (state variable: m(-1) in FOCs)
P_e         % composite energy price (retailer Lagrange multiplier)
P_c         % clean energy price = w_u/A_c
P_d         % dirty energy price = er*P_oil + tau
er          % real exchange rate

% Production (10)
Y_g         % gross output (before damages)
Y_net       % output net de damages = (1-Dam)*Y_g
r_k         % capital rental rate
L_a         % CES labor aggregate
L_s_y       % skilled labor in final production
L_u_y       % unskilled labor in final production
L_u_c       % unskilled labor in clean energy (green jobs)
K_c         % capital in the clean-energy sector (V30: K = K_y + K_c)
L_R_c       % skilled labor in green R&D
w_s         % skilled wage
w_u         % unskilled wage
W_a         % average wage of the CES aggregate
V_KL        % capital-labor composite of final production (nested CES)

% Innovation & growth (2)
A_c         % green technology stock (stationary, SS: delta_c*A_c_ss = eta_c*L_R_c_ss)
gamma_agg   % taux de growth semi-endogene (reporting : gz*Y_net/Y_net(-1))

% Climate (3)
e_em        % domestic emissions = phi_d*E_d
X           % atmospheric carbon stock
Dam         % climate damages Dam = d0+d1*X+d2*X^2

% Carbon-tax recycling (4)
R_tau       % carbon-tax revenue = tau*E_d
T_H         % transfers to HtM
T_R         % transfers to Ricardians
C_n         % non-energy consumption aggregate = lambda*C_H_n + (1-lambda)*C_R_n

% Foreign trade (1)
EX          % exports (pay the oil bill)

% Welfare (recursions, for CEV) (2)
W_H         % HtM welfare (discounted sum)
W_R         % Ricardian welfare
;

%%
%%%%%%%%%%%%% Exogenous variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo
tau         % carbon tax (deterministic Net-Zero trajectory)
P_oil       % world oil price (AR1 or trajectory)
e_ROW       % rest-of-world emissions (exogenous)
theta_policy % R_tau share to R&D; benchmark at t=0, reform from t=1
phi_policy   % HtM share of transfer envelope; benchmark at t=0, reform from t=1
;

%%
%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters
% Preferences
beta        % discount factor
sigma       % CRRA (EIS = 1/sigma)
psi         % inverse Frisch elasticity
chi         % labor weight (calibrated: L_ss = 1/3)
lambda_htM  % HtM share (Campbell-Mankiw)
omega       % energy weight in CES (energy budget share)
epsilon_E   % energy/non-energy substitution elasticity (< 1: rigid energy)
ebar        % energy subsistence floor (Stone-Geary)

% KLE production (shares sum to 1)
alpha_y     % capital share (= 0.30)
beta_y      % labor share (= 0.63)
% energy share = 1-alpha_y-beta_y = 0.07 (implicit)
alpha_c     % capital share in clean energy (V30: Cobb-Douglas K_c^alpha_c * L_u_c^(1-alpha_c))
sigma_L     % skilled/unskilled substitution elasticity (Katz-Murphy)
omega_L     % skilled weight in CES labor (calibrated: w_s/w_u = 1.60)
delta       % capital depreciation rate
phi_k       % investment adjustment cost (I/K ratio)
ass         % neutral TFP (calibrated for Y_ss = 1)

% Energy / mix
alpha_E     % clean-energy weight in energy CES
rho_E       % elasticite substitution clean/dirty
phi_mix     % mix-adjustment friction (greenflation)

% Innovation & growth
gz          % BGP growth rate (world frontier, = 1.003)
delta_c     % green-technology depreciation rate
eta_c       % green R&D productivity (calibrated: delta_c*A_c_ss = eta_c*L_R_c_ss)

% Climate
phi_d       % emissions per unit of dirty energy
delta_clim  % atmospheric carbon absorption
d0 d1 d2   % damage-function parameters (quadratic DICE)

% Open economy
eta_X       % export price elasticity (Henriet et al.)
D_ex_ss     % foreign demand (calibrated for er_ss = 1)

% Steady state (for initval + diagnostics)
gz_ss       % = gz (identity)
tau_ss      % initial carbon tax (= tau_init)
P_oil_ss    % SS oil price
A_c_ss      % SS green technology

% Final-production CES (Schubert specification)
alpha_KL    % capital share in the inner CD composite V = Kprod^alpha_KL * L_a^(1-alpha_KL)
m_Y         % weight of KL composite in outer CES (calibrated at SS)
sigma_Y     % (KL)/energy substitution elasticity in the outer CES
;

%%
%%%%%%%%%%%%% MATLAB SECTION: calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% STRUCTURAL PARAMETERS (directly from the V7 table)
% =========================================================================
issue1_root = getenv('ETANK_ISSUE1_ROOT');
if isempty(issue1_root); issue1_root = pwd; end
issue1_results_dir = fullfile(issue1_root, 'results_v19_issue5');
if ~exist(issue1_results_dir, 'dir'); mkdir(issue1_results_dir); end
issue1_benchmark_file = getenv('ETANK_ISSUE5_BENCHMARK_FILE');
if isempty(issue1_benchmark_file)
    issue1_benchmark_file = fullfile(issue1_results_dir, 'benchmark_calibration.mat');
end

if exist(issue1_benchmark_file, 'file')
    load(issue1_benchmark_file, 'Calib', 'SS0');
    issue1_calib_fields = fieldnames(Calib);
    for issue1_i = 1:numel(issue1_calib_fields)
        issue1_name = issue1_calib_fields{issue1_i};
        eval([issue1_name ' = Calib.(issue1_name);']);
    end
    issue1_ss_fields = fieldnames(SS0);
    for issue1_i = 1:numel(issue1_ss_fields)
        issue1_name = issue1_ss_fields{issue1_i};
        if isnumeric(SS0.(issue1_name)) && isscalar(SS0.(issue1_name))
            eval([issue1_name '_ss = SS0.(issue1_name);']);
        end
    end
    K_total_ss = SS0.K;
    fprintf('\n=== ISSUE 5: loaded selected frozen economy calibration ===\n');
else
fprintf('\n=== ISSUE 3: creating provisional benchmark calibration ===\n');
beta       = 0.997;      % discount factor
sigma      = 1.50;       % CRRA (EIS = 0.67)
psi        = 2.0;        % inverse Frisch
lambda_htM = 0.20;       % final benchmark HtM share
omega      = 0.012;      % fixed benchmark preference parameter (not calibrated)
epsilon_E  = 0.30;       % energy/non-energy substitution (< 1)
alpha_y    = 0.30;       % capital share
beta_y     = 0.63;       % labour share
% energy share = 1 - 0.30 - 0.63 = 0.07 (check: sum = 1)
alpha_c    = 0.20;       % V30: capital share in clean energy (Cobb-Douglas)
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

% --- Final-production CES (Schubert specification) ---
sigma_Y    = 0.50;                        % benchmark: (KL)/energy complementarity
alpha_KL   = alpha_y/(alpha_y+beta_y);   % = 0.3226 (preserves internal K/L split)
% m_Y calibrated in SS loop via energy FOC (Step C); sigma_Y->1 => m_Y -> 0.93
eta_X      = 0.60;       % export price elasticity (Henriet et al.)
theta_benchmark = 0.20;  % benchmark: 20% R_tau to R&D
phi_benchmark   = 0.80;  % benchmark: 80% transfers to HtM

% Policy : initial and terminal carbon tax
tau_ss     = 0.05;       % carbon tax at the initial SS (tau_init)
P_oil_ss   = 1.0;        % normalized oil price

% =========================================================================
% ITERATIVE SS CALIBRATION
% =========================================================================
fprintf('\n=== CALIBRATION V30-CES ===\n');

% --- Step 1 : capital rate of return ---
r_k_ss   = gz^sigma / beta - (1-delta);
fprintf('r_k_ss     = %.6f  (target: real r ~1.5%/quarter)\n', r_k_ss);

% --- Step 2 : output normalization ---
Y_net_ss = 1.0;   % normalisation

% --- Step 3 : capital factor and investment ---
K_ss     = alpha_y * Y_net_ss / r_k_ss;
I_ss     = (gz - (1-delta)) * K_ss;     % from the accumulation equation: gz*K = (1-delta)*K + I
fprintf('K_ss       = %.6f\n', K_ss);
fprintf('I_ss       = %.6f  (I/Y = %.2f%%)\n', I_ss, 100*I_ss/Y_net_ss);

% --- Step 4: energy in production ---
% FOC energy : P_e*E_y = (1-alpha_y-beta_y)*Y_net
energy_share = 1 - alpha_y - beta_y;   % = 0.07
fprintf('energy_share = %.2f (= 1 - %.2f - %.2f)  sum = 1\n', energy_share, alpha_y, beta_y);

% Need the SS energy price to calibrate E_y_ss
% Clean energy price : P_c = w_u/A_c, dirty price : P_d = P_oil + tau_ss
% At SS without mix friction (m constant), the CES gives P_e from P_c and P_d
% First: exchange rate er_ss (calibrate D_ex_ss after), set er_ss = 1
er_ss    = 1.0;
P_d_ss   = er_ss * P_oil_ss + tau_ss;   % = 1.05

% Initialize an iterative loop for the labor and energy markets
A_c_ss   = 1.0;   % green-technology normalization (target)

% Initial guess
L_H_u_ss = 1/3;
L_R_s_ss = 1/3;

% Initial guess for convergence variables (BEFORE the loop)
R_tau_ss  = tau_ss * 0.004;
    L_R_c_ss  = theta_benchmark * R_tau_ss / 1.5;
L_u_c_ss  = 0.005;
    T_H_ss    = phi_benchmark*(1-theta_benchmark)*R_tau_ss/lambda_htM;
    T_R_ss    = (1-phi_benchmark)*(1-theta_benchmark)*R_tau_ss/(1-lambda_htM);

% W_a_ss initialization: CD approximation for first CES loop iteration
W_a_ss = beta_y * Y_net_ss / 0.47;   % L_a ~ 0.47 (initial guess)
m_Y    = alpha_y + beta_y;            % = 0.93 (initial CD value)

% --- Pre-calibration fonction damages (target fixe 1% GDP a +2degC, X=790 GtCO2) ---
% CRITICAL NOTE (CES-spec bug correction): the specification document assumes Dam_ss=0 at the initial SS
% but the actual model has Dam_ss ~ 0.133 (X_ss~2881 GtCO2, pre-existing
% atmospheric carbon stock already accumulated). One must therefore use Y_g_ss =
% Y_net_ss/(1-Dam_ss) -- not Y_net_ss -- in the CES calibration FOCs
% (Steps C/D below). Dam_ss depends on E_d_ss which is only known at the end
% of THIS iteration: we therefore use Dam_ss_prev (lagged one iteration),
% updated at the end of the loop once E_d_ss is computed.
X_2deg   = 790;
d2       = 0.01 / X_2deg^2;
d1       = 0;
d0       = 0;
e_ROW_ss = 10;
Dam_ss_prev = 0.0;   % initial guess (converges to Dam_ss over the iterations)
rho_Y = 1 - 1/sigma_Y;   % CES parameter (Y_g = ass*(m_Y*V^rho+(1-m_Y)*E^rho)^(1/rho))
ass_prev = 1.0;   % initial guess for the TFP factor ass^rho in the FOCs (converges within loop)

for iter_main = 1:500
    % Total labor in production
    L_s_y_ss  = (1-lambda_htM)*L_R_s_ss - L_R_c_ss;
    L_u_y_ss  = lambda_htM*L_H_u_ss - L_u_c_ss;
    if L_s_y_ss <= 0; L_s_y_ss = 0.01; end
    if L_u_y_ss <= 0; L_u_y_ss = 0.01; end

    % CES labor aggregate (CES sigma_L)
    % L_a = [omega_L*L_s_y^((sigma_L-1)/sigma_L) + (1-omega_L)*L_u_y^((sigma_L-1)/sigma_L)]^(sigma_L/(sigma_L-1))
    % To calibrate omega_L: target skill premium w_s/w_u = 1.60
    % FOC : w_s/w_u = (omega_L/(1-omega_L)) * (L_u_y/L_s_y)^(1/sigma_L)
    % => omega_L = (w_s/w_u) * (L_s_y/L_u_y)^(1/sigma_L) / (1 + (w_s/w_u)*(L_s_y/L_u_y)^(1/sigma_L))
    skill_premium_target = 1.60;
    ratio_Lu_Ls = L_u_y_ss / L_s_y_ss;
    omega_L    = skill_premium_target * ratio_Lu_Ls^(-1/sigma_L) / (1 + skill_premium_target * ratio_Lu_Ls^(-1/sigma_L));
    nu_L = (sigma_L-1)/sigma_L;
    L_a_ss = (omega_L*L_s_y_ss^nu_L + (1-omega_L)*L_u_y_ss^nu_L)^(1/nu_L);

    % Wages (W_a_ss from the previous iteration -- CES update downstream)
    w_s_ss = W_a_ss * omega_L * (L_a_ss/L_s_y_ss)^(1/sigma_L);
    w_u_ss = W_a_ss * (1-omega_L) * (L_a_ss/L_u_y_ss)^(1/sigma_L);

    % Clean energy price (V30: Cobb-Douglas unit cost)
    P_c_ss = (1/A_c_ss) * (r_k_ss/alpha_c)^alpha_c * (w_u_ss/(1-alpha_c))^(1-alpha_c);
    % CES energy price (retailer SS without friction)
    % At SS: m constant => no friction term
    % P_e FOC clean  : P_e*alpha_E*(E/E_c)^(1/rho) = P_c
    % P_e FOC dirty    : P_e*(1-alpha_E)*(E/E_d)^(1/rho) = P_d
    % CES minimum-price formula :
    P_e_ss  = (alpha_E^rho_E * P_c_ss^(1-rho_E) + (1-alpha_E)^rho_E * P_d_ss^(1-rho_E))^(1/(1-rho_E));

    % Energy in production (target share 7%, preserved)
    E_y_ss  = energy_share * Y_net_ss / P_e_ss;

    % --- Final-production CES, Steps A-D ---
    % Step A: KL composite
    V_KL_ss = K_ss^alpha_KL * L_a_ss^(1-alpha_KL);
    % Y_g_ss approx (uses Dam_ss lagged from the previous iteration;
    % converges to the true value since Dam_ss is updated each round)
    Y_g_ss_approx = Y_net_ss / max(0.001, 1 - Dam_ss_prev);
    % Step C: m_Y via the energy FOC (CES) -- standard MPK/MPE: dY_g/dX = a*ass^rho*(Y_g/X)^(1/sigma_Y)
    % => P_e = (1-m_Y)*ass_prev^rho_Y*(Y_g/E_y)^(1/sigma_Y)*(1-Dam)  (ass lagged, converges within loop)
    m_Y = 1 - P_e_ss / ( ass_prev^rho_Y * (Y_g_ss_approx / E_y_ss)^(1/sigma_Y) * max(0.001, 1 - Dam_ss_prev) );
    m_Y = max(0.01, min(0.9999, m_Y));
    % Inner K_ss: CES capital-FOC fixed point (r_k = alpha_KL*m_Y*ass^rho*(Y_g/V)^(1/sigma_Y)*(V/K)*(1-Dam))
    for iter_K = 1:100
        V_tmp = K_ss^alpha_KL * L_a_ss^(1-alpha_KL);
        K_new = alpha_KL * m_Y * ass_prev^rho_Y * (Y_g_ss_approx/V_tmp)^(1/sigma_Y) * V_tmp * (1-Dam_ss_prev) / r_k_ss;
        err_K = abs(K_new - K_ss);
        K_ss  = 0.5*K_ss + 0.5*K_new;
        if err_K < 1e-12; break; end
    end
    V_KL_ss = K_ss^alpha_KL * L_a_ss^(1-alpha_KL);
    I_ss    = (gz - (1-delta)) * K_ss;   % updated with CES K_y
    % Step D: aggregate CES wage (W_a = (1-alpha_KL)*m_Y*ass^rho*(Y_g/V)^(1/sigma_Y)*(V/L_a)*(1-Dam))
    W_a_ss_new  = (1-alpha_KL) * m_Y * ass_prev^rho_Y * (Y_g_ss_approx/V_KL_ss)^(1/sigma_Y) * V_KL_ss * (1-Dam_ss_prev) / L_a_ss;
    W_a_ss_prev = W_a_ss;
    W_a_ss      = 0.8*W_a_ss + 0.2*W_a_ss_new;   % damped update

    % Clean and dirty energy (CES allocation)
    % E_c/E_d = (alpha_E/(1-alpha_E))^rho_E * (P_d/P_c)^rho_E
    ratio_Ec_Ed = (alpha_E/(1-alpha_E))^rho_E * (P_d_ss/P_c_ss)^rho_E;
    % E_c = ratio_Ec_Ed * E_d
    % E = [alpha_E*E_c^nu_rho + (1-alpha_E)*E_d^nu_rho]^(1/nu_rho) where nu_rho = (rho-1)/rho
    % => E_d from E_y share (SS simplification : E ~= E_y + E_households)
    % Household energy demand (Stone-Geary CES) -- bootstrap
    P_Z_ss  = ((1-omega) + omega*P_e_ss^(1-epsilon_E))^(1/(1-epsilon_E));

    % Transfers (using R_tau_ss from the previous iteration)
    T_H_ss     = phi_benchmark * (1-theta_benchmark) * R_tau_ss / lambda_htM;
    T_R_ss     = (1-phi_benchmark) * (1-theta_benchmark) * R_tau_ss / (1-lambda_htM);
    L_R_c_ss   = theta_benchmark * R_tau_ss / w_s_ss;

    % HtM budget -- target L_H_u = 1/3 (chi calibrated after the loop)
    ebar = 0.044;  % V30: energy subsistence floor (vs 0.05 in V29)
    L_H_u_ss = 1/3;
    budget_H  = w_u_ss*L_H_u_ss + T_H_ss - P_e_ss*ebar;
    if budget_H <= 0; ebar = 0.5*(w_u_ss*L_H_u_ss + T_H_ss)/P_e_ss; budget_H = w_u_ss*L_H_u_ss + T_H_ss - P_e_ss*ebar; end
    Z_H_ss = budget_H / P_Z_ss;

    % Ricardian budget -- target L_R_s = 1/3
    L_R_s_ss = 1/3;
    income_R  = w_s_ss*L_R_s_ss + r_k_ss*K_ss + T_R_ss - P_e_ss*ebar - I_ss;
    if income_R <= 0; income_R = 0.1; end
    Z_R_ss    = income_R / P_Z_ss;
    L_R_s_new = 1/3;

    % Household energy (Stone-Geary demand)
    C_H_e_ss = ebar + omega*(P_Z_ss/P_e_ss)^epsilon_E * Z_H_ss;
    C_R_e_ss = ebar + omega*(P_Z_ss/P_e_ss)^epsilon_E * Z_R_ss;
    C_H_n_ss = (1-omega)*P_Z_ss^epsilon_E * Z_H_ss;
    C_R_n_ss = (1-omega)*P_Z_ss^epsilon_E * Z_R_ss;

    % Total energy
    E_hh_ss = lambda_htM*C_H_e_ss + (1-lambda_htM)*C_R_e_ss;
    E_ss    = E_hh_ss + E_y_ss;
    % Decomposition into E_c, E_d
    E_d_ss  = E_ss / (1 + ratio_Ec_Ed);
    E_c_ss  = ratio_Ec_Ed * E_d_ss;
    m_ss    = E_c_ss / E_d_ss;   % = ratio_Ec_Ed

    % Update Dam_ss_prev (for the CES block of the next iteration):
    % carbon balance consistent with the current E_d_ss
    X_ss_iter   = (phi_d * E_d_ss + e_ROW_ss) / delta_clim;
    Dam_ss_prev = d0 + d1*X_ss_iter + d2*X_ss_iter^2;

    % Update ass_prev (TFP^rho factor in the CES FOCs, damped fixed point):
    % ass = Y_g_ss / (m_Y*V_KL^rho + (1-m_Y)*E_y^rho)^(1/rho) (definition residuelle CES)
    Y_g_ss_iter = Y_net_ss / max(0.001, 1 - Dam_ss_prev);
    S_iter      = m_Y*V_KL_ss^rho_Y + (1-m_Y)*E_y_ss^rho_Y;
    ass_new     = Y_g_ss_iter / S_iter^(1/rho_Y);
    ass_prev    = 0.5*ass_prev + 0.5*ass_new;

    % Update L_u_c (V30: labor factor demand, not a direct production function)
    L_u_c_ss_new = (1-alpha_c) * P_c_ss * E_c_ss / w_u_ss;

    % Tax revenue and allocations
    R_tau_ss_new = tau_ss * E_d_ss;
    T_H_ss_new   = phi_benchmark*(1-theta_benchmark)*R_tau_ss_new / lambda_htM;
    T_R_ss_new   = (1-phi_benchmark)*(1-theta_benchmark)*R_tau_ss_new / (1-lambda_htM);
    L_R_c_ss_new = theta_benchmark * R_tau_ss_new / w_s_ss;

    % Convergence
    err_main = max(abs([L_u_c_ss_new - L_u_c_ss, R_tau_ss_new - R_tau_ss, ...
                        L_R_c_ss_new - L_R_c_ss, W_a_ss_new - W_a_ss_prev]));

    % L_R_s_ss stays fixed at 1/3 (target)
    L_u_c_ss = 0.5*L_u_c_ss + 0.5*L_u_c_ss_new;
    R_tau_ss = 0.5*R_tau_ss + 0.5*R_tau_ss_new;
    L_R_c_ss = 0.5*L_R_c_ss + 0.5*L_R_c_ss_new;
    T_H_ss   = T_H_ss_new;
    T_R_ss   = T_R_ss_new;

    if err_main < 1e-9; fprintf('SS convergence reached in %d iterations\n', iter_main); break; end
end

% --- V30: Capital in clean energy and total ---
K_c_ss     = alpha_c * P_c_ss * E_c_ss / r_k_ss;
K_total_ss = K_ss + K_c_ss;            % K_ss = K_y (final good only)
I_ss       = (gz - (1-delta)) * K_total_ss;   % investment couvre K_y + K_c
% Ricardian budget corrected (capital income = r_k * K_total)
income_R_corr = w_s_ss*L_R_s_ss + r_k_ss*K_total_ss + T_R_ss - P_e_ss*ebar - I_ss;
if income_R_corr > 0
    Z_R_ss = income_R_corr / P_Z_ss;
end
fprintf('K_c_ss     = %.6f (V30: capital clean energy)\n', K_c_ss);
fprintf('K_total_ss = %.6f (K_y + K_c)\n', K_total_ss);
fprintf('I_ss (V30) = %.6f (couvre K_y + K_c)\n', I_ss);

% --- Calibration of chi (labor weight, target L_ss = 1/3 for each type) ---
% HtM labor FOC : chi*L_H_u^psi = Z_H^(-sigma)*w_u/P_Z
chi_H = Z_H_ss^(-sigma)*w_u_ss / (P_Z_ss * L_H_u_ss^psi);
% Ricardian labor FOC : chi*L_R_s^psi = Z_R^(-sigma)*w_s/P_Z
chi_R = Z_R_ss^(-sigma)*w_s_ss / (P_Z_ss * L_R_s_ss^psi);
chi   = 0.5*(chi_H + chi_R);  % average (should be close if well calibrated)
fprintf('chi (H)    = %.6f | chi (R) = %.6f | chi retenu = %.6f\n', chi_H, chi_R, chi);

% --- Calibration de ebar (target burden HtM ~19%) ---
% burden_H = P_e*C_H_e / (C_H_n + P_e*C_H_e)
% ebar is adjusted until burden_H ~ 0.19 is reached
burden_H_actual = P_e_ss*C_H_e_ss / (C_H_n_ss + P_e_ss*C_H_e_ss);
burden_R_actual = P_e_ss*C_R_e_ss / (C_R_n_ss + P_e_ss*C_R_e_ss);
regressivity    = burden_H_actual / burden_R_actual;
fprintf('burden_H   = %.4f (target ~0.19)\n', burden_H_actual);
fprintf('burden_R   = %.4f (target ~0.07)\n', burden_R_actual);
fprintf('regressivity = %.4f (target ~2.5)\n', regressivity);
fprintf('[Note: adjust ebar manually if burden_H is too far from 0.19]\n');

% --- Calibration of eta_c (SS green innovation: delta_c*A_c_ss = eta_c*L_R_c_ss) ---
eta_c = delta_c * A_c_ss / L_R_c_ss;
fprintf('eta_c      = %.6f (calibrated EE : delta_c*A_c_ss/L_R_c_ss)\n', eta_c);

% --- Calibration of the damage function (MUST be done BEFORE ass) ---
% Target: 1% GDP at +2 deg C (X = 790 GtCO2)
X_2deg = 790;
d2 = 0.01 / X_2deg^2;   % approx 1.603e-8
d1 = 0;
d0 = 0;
fprintf('Dam a +2degC = %.4f (target 1%% GDP, X=790 GtCO2)\n', d0+d1*X_2deg+d2*X_2deg^2);

% --- Calibration of ass (neutral TFP) ---
% X_ss consistent with the model's carbon balance: delta_clim*X = e_em + e_ROW
e_ROW_ss = 10;
X_ss = (phi_d * E_d_ss + e_ROW_ss) / delta_clim;  % ~2891 GtCO2
Dam_ss = d0 + d1*X_ss + d2*X_ss^2;
Y_g_ss = Y_net_ss / (1 - Dam_ss);
% Step E: calibrate ass (neutral TFP, CES formula)
V_KL_ss = K_ss^alpha_KL * L_a_ss^(1-alpha_KL);   % final value after convergence
ass = Y_g_ss / ( m_Y*V_KL_ss^(1-1/sigma_Y) + (1-m_Y)*E_y_ss^(1-1/sigma_Y) )^(1/(1-1/sigma_Y));
fprintf('ass        = %.6f (normalized CES TFP)\n', ass);
% Step F : ex post verification of the capital FOC (deviation < 1e-6 required)
r_k_check = alpha_KL * m_Y * ass^(1-1/sigma_Y) * (Y_g_ss/V_KL_ss)^(1/sigma_Y) * V_KL_ss*(1-Dam_ss)/K_ss;
fprintf('r_k FOC check : r_k_ss=%.6f vs CES=%.6f (ecart=%.2e)\n', r_k_ss, r_k_check, abs(r_k_ss-r_k_check));
fprintf('m_Y        = %.6f | alpha_KL = %.6f | sigma_Y = %.4f\n', m_Y, alpha_KL, sigma_Y);
fprintf('V_KL_ss    = %.6f | K_ss CES = %.6f (vs K_ss CD = %.6f)\n', V_KL_ss, K_ss, alpha_y*Y_net_ss/r_k_ss);
fprintf('Dam_ss     = %.6f (Dam SS, X_ss=%.1f GtCO2)\n', Dam_ss, X_ss);
fprintf('X_ss       = %.1f (SS consistent with the carbon balance)\n', X_ss);

% --- Calibration of D_ex_ss (foreign demand, target er_ss=1) ---
D_ex_ss = P_oil_ss * er_ss * E_d_ss * er_ss^eta_X;
fprintf('D_ex_ss    = %.6f (SS exports = oil bill)\n', D_ex_ss);

% --- Welfare SS (for CEV) ---
% U_H_ss = Z_H_ss^(1-sigma)/(1-sigma) - chi*L_H_u_ss^(1+psi)/(1+psi)
U_H_ss = Z_H_ss^(1-sigma)/(1-sigma) - chi*L_H_u_ss^(1+psi)/(1+psi);
U_R_ss = Z_R_ss^(1-sigma)/(1-sigma) - chi*L_R_s_ss^(1+psi)/(1+psi);
% W_ss = U_ss / (1 - beta*gz^(1-sigma))
disc_W = beta * gz^(1-sigma);
if abs(1-disc_W) < 1e-8; error('Welfare discount factor = 1, instable'); end
W_H_ss = U_H_ss / (1 - disc_W);
W_R_ss = U_R_ss / (1 - disc_W);
fprintf('W_H_ss     = %.4f | W_R_ss = %.4f\n', W_H_ss, W_R_ss);
fprintf('Welfare discount factor : beta*gz^(1-sigma) = %.6f\n', disc_W);

% --- Real rate SS (for initval) ---
rr_ss = gz^sigma / beta;  % short-term gross real rate
r_ss  = rr_ss;            % real model (no inflation)

% --- SS aggregates ---
C_n_ss = lambda_htM*C_H_n_ss + (1-lambda_htM)*C_R_n_ss;
EX_ss  = P_oil_ss * er_ss * E_d_ss;
Lambda_R_ss = Z_R_ss^(-sigma) / P_Z_ss;

% --- Tobin's q SS ---
q_K_ss = 1 / (1 - phi_k*(I_ss/K_total_ss - delta));   % V30 : K_total
gamma_agg_ss = gz;

fprintf('\n--- SS SUMMARY V30-CES ---\n');
fprintf('Y_net_ss   = %.4f | K_ss     = %.4f | I_ss    = %.4f\n', Y_net_ss, K_ss, I_ss);
fprintf('w_u_ss     = %.4f | w_s_ss   = %.4f | prem.sal = %.4f\n', w_u_ss, w_s_ss, w_s_ss/w_u_ss);
fprintf('L_H_u_ss   = %.4f | L_R_s_ss = %.4f\n', L_H_u_ss, L_R_s_ss);
fprintf('P_e_ss     = %.4f | P_c_ss   = %.4f | P_d_ss  = %.4f\n', P_e_ss, P_c_ss, P_d_ss);
fprintf('E_ss       = %.6f | E_c_ss   = %.6f | E_d_ss  = %.6f\n', E_ss, E_c_ss, E_d_ss);
fprintf('A_c_ss     = %.4f | L_R_c_ss = %.6f | L_u_c_ss= %.6f\n', A_c_ss, L_R_c_ss, L_u_c_ss);
fprintf('R_tau_ss   = %.6f | T_H_ss   = %.6f | T_R_ss  = %.6f\n', R_tau_ss, T_H_ss, T_R_ss);
fprintf('er_ss      = %.4f | EX_ss    = %.6f\n', er_ss, EX_ss);
fprintf('P_Z_ss     = %.4f | Z_H_ss   = %.4f | Z_R_ss  = %.4f\n', P_Z_ss, Z_H_ss, Z_R_ss);
fprintf('omega_L    = %.4f (calibrated skill premium)\n', omega_L);
fprintf('chi        = %.4f (calibrated L_ss~1/3)\n', chi);
fprintf('eta_c      = %.6f (calibrated A_c_ss=1)\n', eta_c);
fprintf('D_ex_ss    = %.6f (calibrated er_ss=1)\n', D_ex_ss);
fprintf('=====================\n');

% Freeze the benchmark as the single source of truth for every scenario.
Calib = struct();
Calib.beta=beta; Calib.sigma=sigma; Calib.psi=psi; Calib.chi=chi;
Calib.lambda_htM=lambda_htM; Calib.omega=omega; Calib.epsilon_E=epsilon_E;
Calib.ebar=ebar; Calib.alpha_y=alpha_y; Calib.beta_y=beta_y;
Calib.alpha_c=alpha_c; Calib.sigma_L=sigma_L; Calib.omega_L=omega_L;
Calib.delta=delta; Calib.phi_k=phi_k; Calib.ass=ass; Calib.A=ass;
Calib.alpha_E=alpha_E; Calib.rho_E=rho_E; Calib.phi_mix=phi_mix;
Calib.gz=gz; Calib.delta_c=delta_c; Calib.eta_c=eta_c;
Calib.phi_d=phi_d; Calib.delta_clim=delta_clim;
Calib.d0=d0; Calib.d1=d1; Calib.d2=d2;
Calib.eta_X=eta_X; Calib.D_ex_ss=D_ex_ss;
Calib.gz_ss=gz_ss; Calib.tau_ss=tau_ss; Calib.P_oil_ss=P_oil_ss;
Calib.A_c_ss=A_c_ss; Calib.alpha_KL=alpha_KL; Calib.m_Y=m_Y;
Calib.sigma_Y=sigma_Y; Calib.theta0=theta_benchmark;
Calib.phi0=phi_benchmark; Calib.tau0=tau_ss;

SS0 = struct();
SS0.P_Z=P_Z_ss; SS0.Z_H=Z_H_ss; SS0.Z_R=Z_R_ss;
SS0.C_H_n=C_H_n_ss; SS0.C_R_n=C_R_n_ss;
SS0.C_H_e=C_H_e_ss; SS0.C_R_e=C_R_e_ss;
SS0.L_H_u=L_H_u_ss; SS0.L_R_s=L_R_s_ss;
SS0.Lambda_R=Lambda_R_ss; SS0.I=I_ss; SS0.K=K_total_ss;
SS0.q_K=q_K_ss; SS0.E=E_ss; SS0.E_c=E_c_ss; SS0.E_d=E_d_ss;
SS0.E_y=E_y_ss; SS0.m=m_ss; SS0.P_e=P_e_ss; SS0.P_c=P_c_ss;
SS0.P_d=P_d_ss; SS0.er=er_ss; SS0.Y_g=Y_g_ss; SS0.Y_net=Y_net_ss;
SS0.r_k=r_k_ss; SS0.L_a=L_a_ss; SS0.L_s_y=L_s_y_ss;
SS0.L_u_y=L_u_y_ss; SS0.L_u_c=L_u_c_ss; SS0.K_c=K_c_ss;
SS0.L_R_c=L_R_c_ss; SS0.w_s=w_s_ss; SS0.w_u=w_u_ss;
SS0.W_a=W_a_ss; SS0.V_KL=V_KL_ss; SS0.A_c=A_c_ss;
SS0.gamma_agg=gamma_agg_ss; SS0.e_em=phi_d*E_d_ss; SS0.X=X_ss;
SS0.Dam=Dam_ss; SS0.R_tau=R_tau_ss; SS0.T_H=T_H_ss; SS0.T_R=T_R_ss;
SS0.C_n=C_n_ss; SS0.EX=EX_ss; SS0.W_H=W_H_ss; SS0.W_R=W_R_ss;
SS0.tau=tau_ss; SS0.P_oil=P_oil_ss; SS0.e_ROW=e_ROW_ss;
SS0.theta_policy=theta_benchmark; SS0.phi_policy=phi_benchmark;
SS0.burden_H=burden_H_actual; SS0.burden_R=burden_R_actual;
SS0.burden_ratio=regressivity;
save(issue1_benchmark_file, 'Calib', 'SS0');
fprintf('Saved unique benchmark: %s\n', issue1_benchmark_file);
end

% Loading MATLAB variables is not sufficient for Dynare: synchronize every
% frozen structural parameter explicitly with M_.params before steady/solver.
issue1_param_names = cellstr(M_.param_names);
for issue1_i = 1:numel(issue1_param_names)
    issue1_name = strtrim(issue1_param_names{issue1_i});
    if isfield(Calib, issue1_name)
        set_param_value(issue1_name, Calib.(issue1_name));
    end
end

%%
%%%%%%%%%%%%% MODEL BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model;
% =========================================================================
% BLOCK 1: STONE-GEARY PREFERENCES (CES)
% =========================================================================

% [1] CES composite price index, energy/non-energy
% P_Z = [(1-omega) + omega*P_e^(1-epsilon_E)]^(1/(1-epsilon_E))
P_Z = ((1-omega) + omega*P_e^(1-epsilon_E))^(1/(1-epsilon_E));

% [2] HtM non-energy demand (numeraire P_n = 1)
C_H_n = (1-omega)*P_Z^epsilon_E * Z_H;

% [3] HtM energy demand (supernumerary + floor)
% C_H_e = ebar + omega*(P_Z/P_e)^epsilon_E * Z_H
C_H_e = ebar + omega*(P_Z/P_e)^epsilon_E * Z_H;

% [4] Ricardian non-energy demand
C_R_n = (1-omega)*P_Z^epsilon_E * Z_R;

% [5] Ricardian energy demand
C_R_e = ebar + omega*(P_Z/P_e)^epsilon_E * Z_R;

% [6] HtM budget (in terms of Z_H)
% P_Z*Z_H + P_e*ebar = w_u*L_H_u + T_H
% => P_Z*Z_H = w_u*L_H_u + T_H - P_e*ebar
P_Z*Z_H = w_u*L_H_u + T_H - P_e*ebar;

% [7] HtM labor FOC (Stone-Geary: MU of consumption = Z_H^{-sigma}/P_Z)
% chi*L_H_u^psi = Z_H^(-sigma)*w_u/P_Z
chi*L_H_u^psi = Z_H^(-sigma)*w_u/P_Z;

% [8] Ricardian marginal utility of income
Lambda_R = Z_R^(-sigma)/P_Z;

% [9] Ricardian budget -- Issue 2 aggregation correction.
% K and I are economy-wide; each Ricardian owns 1/(1-lambda) of them.
P_Z*Z_R = w_s*L_R_s + (r_k*K(-1)-I)/(1-lambda_htM) + T_R - P_e*ebar;

% [10] Ricardian labor FOC
chi*L_R_s^psi = Z_R^(-sigma)*w_s/P_Z;

% [11] Tobin's q (investment optimality, I/K costs)
% q_K*(1 - phi_k*(I/K(-1) - delta)) = 1
q_K*(1 - phi_k*(I/K(-1) - delta)) = 1;

% [12] Capital accumulation (detrended by gz)
% gz*K = (1-delta)*K(-1) + I - phi_k/2*(I/K(-1) - delta)^2 * K(-1)
gz*K = (1-delta)*K(-1) + I - phi_k/2*(I/K(-1) - delta)^2*K(-1);

% [13] Euler equation (investment + forward-looking Tobin's Q)
% Lambda_R*q_K = beta*gz^(-sigma)*Lambda_R(+1)*(r_k(+1) + q_K(+1)*(1-delta))
Lambda_R*q_K = beta*gz^(-sigma)*Lambda_R(+1)*(r_k(+1) + q_K(+1)*(1-delta));

% =========================================================================
% BLOCK 2: ENERGY RETAILER (CES + mix-adjustment friction)
% =========================================================================

% [14] CES energy aggregator
% E = [alpha_E*E_c^((rho_E-1)/rho_E) + (1-alpha_E)*E_d^((rho_E-1)/rho_E)]^(rho_E/(rho_E-1))
E = (alpha_E*E_c^((rho_E-1)/rho_E) + (1-alpha_E)*E_d^((rho_E-1)/rho_E))^(rho_E/(rho_E-1));

% [15] Energy mix ratio (state variable: m(-1) in the FOCs)
m = E_c/E_d;

% [16] Dirty energy price (taxed)
P_d = er*P_oil + tau;

% [17] Clean energy price -- V30: Cobb-Douglas unit cost (vs w_u/A_c in V29)
% P_c = (1/A_c)*(r_k/alpha_c)^alpha_c * (w_u/(1-alpha_c))^(1-alpha_c)
P_c = (1/A_c) * (r_k/alpha_c)^alpha_c * (w_u/(1-alpha_c))^(1-alpha_c);

% [18] Retailer FOC for clean energy (with mix friction)
% P_e*alpha_E*(E/E_c)^(1/rho_E) = P_c + phi_mix*(m - m(-1))*P_e*E/E_d
P_e*alpha_E*(E/E_c)^(1/rho_E) = P_c + phi_mix*(m - m(-1))*P_e*E/E_d;

% [19] Retailer FOC for dirty energy (with mix friction)
% P_e*(1-alpha_E)*(E/E_d)^(1/rho_E) = P_d - phi_mix*(m - m(-1))*E_c*P_e*E/E_d^2
P_e*(1-alpha_E)*(E/E_d)^(1/rho_E) = P_d - phi_mix*(m - m(-1))*E_c*P_e*E/E_d^2;

% =========================================================================
% BLOCK 3: FINAL PRODUCTION (KLE, shares sum to 1)
% =========================================================================

% [20] Capital-labor composite (inner Cobb-Douglas, constant returns)
V_KL = (K(-1)-K_c)^alpha_KL * L_a^(1-alpha_KL);

% [20b] nested CES final production : composite KL + energy
Y_g = ass * ( m_Y*V_KL^(1-1/sigma_Y) + (1-m_Y)*E_y^(1-1/sigma_Y) )^(1/(1-1/sigma_Y));

% [21] Output net of damages (unchanged)
Y_net = (1 - Dam) * Y_g;

% [22] Capital FOC (outer CES, via the V_KL composite)
% NOTE: the factor ass^(1-1/sigma_Y) is required by the standard CES derivative
% MPK = m_Y*ass^rho*(Y_g/V)^(1/sigma_Y); the second term is V_KL/K (not Y_g/K)
r_k = alpha_KL * m_Y * ass^(1-1/sigma_Y) * (Y_g/V_KL)^(1/sigma_Y) * V_KL/(K(-1)-K_c) * (1-Dam);

% [23] Energy FOC in production (outer CES)
P_e = (1-m_Y) * ass^(1-1/sigma_Y) * (Y_g/E_y)^(1/sigma_Y) * (1-Dam);

% [24] Aggregate labor FOC (outer CES, via the V_KL composite)
W_a = (1-alpha_KL) * m_Y * ass^(1-1/sigma_Y) * (Y_g/V_KL)^(1/sigma_Y) * V_KL/L_a * (1-Dam);

% [25] CES labor aggregate (sigma_L)
L_a = (omega_L*L_s_y^((sigma_L-1)/sigma_L) + (1-omega_L)*L_u_y^((sigma_L-1)/sigma_L))^(sigma_L/(sigma_L-1));

% [26] Skilled-labor FOC in production
w_s = W_a*omega_L*(L_a/L_s_y)^(1/sigma_L);

% [27] Unskilled-labor FOC in production
w_u = W_a*(1-omega_L)*(L_a/L_u_y)^(1/sigma_L);

% [28] Unskilled-labor factor demand in clean energy -- V30
% L_u_c = (1-alpha_c)*P_c*E_c/w_u  (labor share of cost, zero profit under CRS)
L_u_c = (1-alpha_c)*P_c*E_c/w_u;

% [46] Capital factor demand in clean energy -- V30 (new equation)
% K_c = alpha_c*P_c*E_c/r_k  (capital share of cost, zero profit under CRS)
K_c = alpha_c*P_c*E_c/r_k;

% =========================================================================
% BLOCK 4: GREEN INNOVATION AND GROWTH
% =========================================================================

% [29] Green technology stock (depreciating knowledge stock)
% A_c = (1-delta_c)*A_c(-1) + eta_c*L_R_c
A_c = (1-delta_c)*A_c(-1) + eta_c*L_R_c;

% [30] Semi-endogenous growth rate (reporting only)
gamma_agg = gz*Y_net/Y_net(-1);

% =========================================================================
% BLOCK 5: CLIMATE (quadratic DICE)
% =========================================================================

% [31] Domestic emissions
e_em = phi_d*E_d;

% [32] Atmospheric carbon stock
X = (1-delta_clim)*X(-1) + e_em + e_ROW;

% [33] DICE damage function
Dam = d0 + d1*X + d2*X^2;

% =========================================================================
% BLOCK 6: GOVERNMENT -- Carbon tax and recycling
% =========================================================================
% ISSUE 1 architecture: theta_policy and phi_policy are deterministic
% exogenous policy variables. Row t=0 is the common benchmark; the reform
% enters in the first simulated period. Structural parameters never move.

% [34] Carbon-tax revenue
R_tau = tau*E_d;

% [35] Green R&D subsidy (L_R_c financing)
L_R_c = theta_policy*R_tau/w_s;

% [36] Transfers to HtM households (per capita)
T_H = phi_policy*(1-theta_policy)*R_tau/lambda_htM;

% [37] Transfers to Ricardians (per capita)
T_R = (1-phi_policy)*(1-theta_policy)*R_tau/(1-lambda_htM);

% =========================================================================
% BLOCK 7: MARKET CLEARING + OPEN ECONOMY
% =========================================================================

% [38] Non-energy consumption aggregate
C_n = lambda_htM*C_H_n + (1-lambda_htM)*C_R_n;

% Note: the goods market (Y_net = C_n + I + adj_cost) is implicit by Walras' law.
% Post-simulation check: residual must be < 1e-6.

% [39] Unskilled labor market
% lambda*L_H_u = L_u_y + L_u_c
lambda_htM*L_H_u = L_u_y + L_u_c;

% [40] Skilled labor market
% (1-lambda)*L_R_s = L_s_y + L_R_c
(1-lambda_htM)*L_R_s = L_s_y + L_R_c;

% [41] Energy market (total demand = household consumption + production)
% E = lambda*C_H_e + (1-lambda)*C_R_e + E_y
E = lambda_htM*C_H_e + (1-lambda_htM)*C_R_e + E_y;

% [42] Foreign export demand
EX = D_ex_ss*er^(-eta_X);

% [43] Open-economy closure (exports = oil bill, trade balance = 0)
% EX = P_oil*er*E_d => imports financed by final-goods exports
EX = P_oil*er*E_d;

% =========================================================================
% BLOCK 8: WELFARE (recursions, for CEV)
% =========================================================================

% [44] HtM welfare (backward recursion of present value)
% W_H = [Z_H^(1-sigma)/(1-sigma) - chi*L_H_u^(1+psi)/(1+psi)] + beta*gz^(1-sigma)*W_H(+1)
W_H = Z_H^(1-sigma)/(1-sigma) - chi*L_H_u^(1+psi)/(1+psi) + beta*gz^(1-sigma)*W_H(+1);

% [45] Ricardian welfare
W_R = Z_R^(1-sigma)/(1-sigma) - chi*L_R_s^(1+psi)/(1+psi) + beta*gz^(1-sigma)*W_R(+1);

end;

%%
%%%%%%%%%%%%% PERFECT FORESIGHT SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Horizon: T_sim = 200 quarters (50 years)
% Carbon tax: convex ramp from tau_init=0.05 to tau_final=0.40 over 120 quarters (30 years)

@#ifndef T_SIM_QUARTERS
@#define T_SIM_QUARTERS = 200
@#endif
@#define T_sim = T_SIM_QUARTERS
@#ifndef TAU_FINAL
@#define TAU_FINAL = 0.40
@#endif
@#ifndef KAPPA_TAX
@#define KAPPA_TAX = 2
@#endif

% Initial conditions (initial SS, tau=0.05)
initval;
tau      = tau_ss;
@#if P_OIL_PATH_MODE == 1
P_oil    = P_oil_ss;
@#else
P_oil    = @{P_OIL_FINAL};
@#endif
e_ROW    = e_ROW_ss;
theta_policy = 0.20;
phi_policy   = 0.80;
er       = er_ss;
P_d      = er_ss*P_oil_ss + tau_ss;
A_c      = A_c_ss;
P_c      = P_c_ss;
P_e      = P_e_ss;
P_Z      = P_Z_ss;
m        = m_ss;
K        = K_total_ss;   % V30 : K_total = K_y + K_c
K_c      = K_c_ss;       % V30 : capital clean energy SS initial
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

% Terminal conditions. TAU_FINAL is the single source used by endval and
% by the full exogenous trajectory below.
endval;
tau      = @{TAU_FINAL};
P_oil    = P_oil_ss;
e_ROW    = e_ROW_ss;
theta_policy = @{THETA_POLICY};
phi_policy   = @{PHI_POLICY};
end;

% Shock: tau_init to tau_final over 120 quarters (convex ramp)
% The exact trajectory is injected further below, directly into this .mod,
% into oo_.exo_simul after perfect_foresight_setup (see the kappa_v=2 block)
shocks;
var tau;
periods 1;
values 0.05;
end;
% Note: the block below (after perfect_foresight_setup) overrides
% oo_.exo_simul with the exact trajectory (convex ramp, kappa=2)

%%
%%%%%%%%%%%%% NUMERICAL SS (after initval/endval/shocks) %%%%%%%%%%%%%%%%%%
% Continue the terminal steady state from tau_0 to TAU_FINAL. This is a
% numerical continuation only; it does not alter equations or tolerances.
homotopy_setup;
tau, tau_ss, @{TAU_FINAL};
@#if P_OIL_PATH_MODE == 0
P_oil, P_oil_ss, @{P_OIL_FINAL};
@#endif
end;
steady(homotopy_mode=1, homotopy_steps=20);
check;

%%
%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perfect_foresight_setup(periods = @{T_sim});

%% Inject the SMOOTH tau trajectory (quadratic ramp) into oo_.exo_simul
% Replaces the abrupt jump from the shocks block (which only set t=1)
tau_idx_v   = find(strcmp(M_.exo_names, 'tau'));
poil_idx_v  = find(strcmp(M_.exo_names, 'P_oil'));
tau_init_v  = 0.05;
tau_final_v = @{TAU_FINAL};
T_carbon_v  = 120;
kappa_v     = @{KAPPA_TAX};
T_sim_v     = @{T_sim};   % = 200
% oo_.exo_simul est (T_sim+2) x n_exo :
%   row 1 = t=0 (initval), rows 2..T+1 = transition, row T+2 = endval
oo_.exo_simul(1, tau_idx_v) = tau_init_v;
oo_.exo_simul(1, poil_idx_v) = P_oil_ss;
theta_idx_v = find(strcmp(M_.exo_names, 'theta_policy'));
phi_idx_v   = find(strcmp(M_.exo_names, 'phi_policy'));
oo_.exo_simul(1, theta_idx_v) = 0.20;
oo_.exo_simul(1, phi_idx_v)   = 0.80;
for t_v = 1:T_sim_v
@#if IDENT_INSTANT_SHOCK == 1
    oo_.exo_simul(t_v+1, tau_idx_v) = tau_final_v;
@#else
    if t_v <= T_carbon_v
        oo_.exo_simul(t_v+1, tau_idx_v) = tau_init_v + ...
            (tau_final_v - tau_init_v) * (t_v / T_carbon_v)^kappa_v;
    else
        oo_.exo_simul(t_v+1, tau_idx_v) = tau_final_v;
    end
@#endif
    oo_.exo_simul(t_v+1, theta_idx_v) = @{THETA_POLICY};
    oo_.exo_simul(t_v+1, phi_idx_v)   = @{PHI_POLICY};
@#if P_OIL_PATH_MODE == 1
    oo_.exo_simul(t_v+1, poil_idx_v)  = P_oil_ss + (@{P_OIL_FINAL}-P_oil_ss)*(@{P_OIL_PERSISTENCE})^(t_v-1);
@#else
    oo_.exo_simul(t_v+1, poil_idx_v)  = @{P_OIL_FINAL};
@#endif
end
@#if P_OIL_PATH_MODE == 1
oo_.exo_simul(T_sim_v+2, poil_idx_v) = P_oil_ss;
@#else
oo_.exo_simul(T_sim_v+2, poil_idx_v) = @{P_OIL_FINAL};
@#endif
oo_.exo_simul(T_sim_v+2, tau_idx_v) = tau_final_v;
oo_.exo_simul(T_sim_v+2, theta_idx_v) = @{THETA_POLICY};
oo_.exo_simul(T_sim_v+2, phi_idx_v)   = @{PHI_POLICY};

perfect_foresight_solver(maxit = 500, tolf = 1e-7,
    homotopy_initial_step_size = 0.001, homotopy_min_step_size = 1e-6);

%%
%%%%%%%%%%%%% RESULTS EXTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    % --- Energy ---
    P_e_sim    = oo_.endo_simul(get_idx('P_e'),    2:T_sim_read+1)';
    E_d_sim    = oo_.endo_simul(get_idx('E_d'),    2:T_sim_read+1)';
    E_c_sim    = oo_.endo_simul(get_idx('E_c'),    2:T_sim_read+1)';
    E_sim      = oo_.endo_simul(get_idx('E'),      2:T_sim_read+1)';
    m_sim      = oo_.endo_simul(get_idx('m'),      2:T_sim_read+1)';

    % --- Distribution and welfare ---
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

    % --- Innovation and green jobs ---
    A_c_sim    = oo_.endo_simul(get_idx('A_c'),      2:T_sim_read+1)';
    L_R_c_sim  = oo_.endo_simul(get_idx('L_R_c'),    2:T_sim_read+1)';
    L_u_c_sim  = oo_.endo_simul(get_idx('L_u_c'),    2:T_sim_read+1)';
    gamma_sim  = oo_.endo_simul(get_idx('gamma_agg'),2:T_sim_read+1)';

    % --- Climate ---
    X_sim      = oo_.endo_simul(get_idx('X'),      2:T_sim_read+1)';
    Dam_sim    = oo_.endo_simul(get_idx('Dam'),    2:T_sim_read+1)';
    e_em_sim   = oo_.endo_simul(get_idx('e_em'),   2:T_sim_read+1)';

    % --- Burdens and regressivity ---
    burden_H_sim = P_e_sim.*C_H_e_sim ./ (C_H_n_sim + P_e_sim.*C_H_e_sim);
    burden_R_sim = P_e_sim.*C_R_e_sim ./ (C_R_n_sim + P_e_sim.*C_R_e_sim);
    regress_sim  = burden_H_sim ./ burden_R_sim;
    prem_sal_sim = w_s_sim ./ w_u_sim;

    % ISSUE 1: normalize against the common pre-announcement benchmark,
    % never against the terminal policy-specific steady state or traj(t=1).
    Y_ss   = SS0.Y_net;
    Pe_ss  = SS0.P_e;
    Ac_ss  = SS0.A_c;
    WH_ss  = SS0.W_H;
    WR_ss  = SS0.W_R;

    fprintf('\n=== V30-CES RESULTS AT 30 YEARS (period 120) ===\n');
    fprintf('Y_net     : %+.2f%% vs SS\n', 100*(Y_net_sim(120)/Y_ss-1));
    fprintf('P_e       : %+.2f%% vs SS\n', 100*(P_e_sim(120)/Pe_ss-1));
    fprintf('A_c       : %+.2f%% vs SS\n', 100*(A_c_sim(120)/Ac_ss-1));
    fprintf('burden_H  : %.4f (SS: %.4f)\n', burden_H_sim(120), P_e_ss*C_H_e_ss/(C_H_n_ss+P_e_ss*C_H_e_ss));
    fprintf('regressiv : %.4f (SS: %.4f)\n', regress_sim(120), SS0.burden_ratio);

    % S-CEV (Supernumerary Consumption Equivalent Variation).
    % Scale the consumption composite only; keep labor disutility fixed.
    % This is the Lucas transform used in the manuscript and grid drivers.
    disc_cev = beta*gz^(1-sigma);
    VC_H_ss = (SS0.Z_H^(1-sigma)/(1-sigma))/(1-disc_cev);
    VL_H_ss = (chi*SS0.L_H_u^(1+psi)/(1+psi))/(1-disc_cev);
    VC_R_ss = (SS0.Z_R^(1-sigma)/(1-sigma))/(1-disc_cev);
    VL_R_ss = (chi*SS0.L_R_s^(1+psi)/(1+psi))/(1-disc_cev);
    WH_t1 = W_H_sim(1);   WR_t1 = W_R_sim(1);
    SCEV_H = 100*(((WH_t1+VL_H_ss)/VC_H_ss)^(1/(1-sigma))-1);
    SCEV_R = 100*(((WR_t1+VL_R_ss)/VC_R_ss)^(1/(1-sigma))-1);
    fprintf('S-CEV HtM : %+.4f%% | S-CEV Ric : %+.4f%%\n', SCEV_H, SCEV_R);

else
    warning('V30-CES simulation did not converge. Check initval/endval and the residuals.');
end
