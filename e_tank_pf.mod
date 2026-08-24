%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% E-TANK PERFECT-FORESIGHT MODEL
%%%% Thesis : Innovate_or_compensate
%%%% Surpervisors: Katheline SCHUBERT, Gauthier VERMANDEL
%%%% Author: Asmae EL MOUHSSINE 
%%%%
%%%% Two-agent environmental general-equilibrium model with:
%%%%   - hand-to-mouth and Ricardian households;
%%%%   - Stone-Geary energy demand;
%%%%   - clean and dirty energy;
%%%%   - endogenous clean-knowledge accumulation;
%%%%   - nested CES final production;
%%%%   - carbon-tax revenue recycling between clean R&D and transfers;
%%%%   - small-open-economy trade closure;
%%%%   - climate damages.
%%%%
%%%% Policy instruments:
%%%%   theta_policy : share of carbon-tax revenue allocated to clean R&D
%%%%   phi_policy   : HtM share of the remaining transfer envelope
%%%%
%%%% Endogenous variables: 47
%%%% Exogenous variables: tau, P_oil, e_ROW, theta_policy, phi_policy
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
%%%%%%%%%%%%% Endogenous variables (47) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var
% Preferences Stone-Geary CES (3)
P_Z         % CES price index energy/non-energy
Z_H         % HtM welfare composite (Stone-Geary CES)
Z_R         % Ricardian welfare composite

% Consumption by type and good (4)
C_H_n       % HtM non-energy consumption
C_R_n       % Ricardian non-energy consumption
C_H_e       % HtM energy consumption (= ebar + supernumerary component)
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
P_c         % clean-energy unit cost
P_d         % dirty energy price = er*P_oil + tau
er          % real exchange rate

% Production (10)
Y_g         % gross output (before damages)
Y_net       % output net of climate damages = (1-Dam)*Y_g
r_k         % capital rental rate
L_a         % CES labor aggregate
L_s_y       % skilled labor in final production
L_u_y       % unskilled labor in final production
L_u_c       % unskilled labor in clean energy (green jobs)
K_c         % capital in the clean-energy sector (K = K_y + K_c)
L_R_c       % skilled labor in green R&D
w_s         % skilled wage
w_u         % unskilled wage
W_a         % average wage of the CES aggregate
V_KL        % capital-labor composite of final production (nested CES)

% Innovation & growth (2)
A_c         % green technology stock (stationary, SS: delta_c*A_c_ss = eta_c*L_R_c_ss)
gamma_agg   % semi-endogenous aggregate growth rate (reported as gz*Y_net/Y_net(-1))

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
tau         % deterministic carbon-tax policy path
P_oil       % world oil price (AR1 or trajectory)
e_ROW       % rest-of-world emissions (normalized model units, exogenous)
theta_policy % R_tau share to R&D; common baseline at t=0, alternative policy from t=1
phi_policy   % HtM transfer share; common baseline at t=0, alternative policy from t=1
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
ebar        % Stone-Geary committed energy component

% Capital-labor-energy production
alpha_y     % capital share (= 0.30)
beta_y      % labor share (= 0.63)
% energy share = 1-alpha_y-beta_y = 0.07 (implicit)
alpha_c     % capital share in clean-energy Cobb-Douglas production
sigma_L     % skilled/unskilled substitution elasticity (Katz-Murphy)
omega_L     % skilled weight in CES labor (calibrated: w_s/w_u = 1.60)
delta       % capital depreciation rate
phi_k       % investment adjustment cost (I/K ratio)
ass         % neutral TFP (calibrated for Y_ss = 1)

% Energy / mix
alpha_E     % clean-energy weight in energy CES
rho_E       % clean/dirty energy substitution elasticity
phi_mix     % clean-dirty energy-mix adjustment friction

% Innovation & growth
gz          % BGP growth rate (world frontier, = 1.003)
delta_c     % green-technology depreciation rate
eta_c       % green R&D productivity (calibrated: delta_c*A_c_ss = eta_c*L_R_c_ss)

% Climate
phi_d       % emissions per unit of dirty energy
delta_clim  % atmospheric carbon absorption
d0 d1 d2   % reduced-form quadratic climate-damage parameters

% Open economy
eta_X       % export price elasticity (Henriet et al.)
D_ex_ss     % foreign demand (calibrated for er_ss = 1)

% Steady state (for initval + diagnostics)
gz_ss       % = gz (identity)
tau_ss      % initial carbon tax (= tau_init)
P_oil_ss    % SS oil price
A_c_ss      % SS green technology

% Final-production nested CES
alpha_KL    % capital share in the inner CD composite V = Kprod^alpha_KL * L_a^(1-alpha_KL)
m_Y         % weight of KL composite in outer CES (calibrated at SS)
sigma_Y     % (KL)/energy substitution elasticity in the outer CES
;

%%
%%%%%%%%%%%%% MATLAB SECTION: calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% FINAL BASELINE CALIBRATION
% =========================================================================
% The final thesis calibration and complete initial steady state are embedded
% below so that the reference experiment is reproducible without an external
% MAT file. An explicit calibration file may still be supplied for robustness
% economies through ETANK_CALIBRATION_FILE. The legacy environment-variable
% name is accepted as a backward-compatible alias for existing drivers.
calibration_file = getenv('ETANK_CALIBRATION_FILE');

if ~isempty(calibration_file)
    if ~exist(calibration_file, 'file')
        error('Specified E-TANK calibration file does not exist: %s', calibration_file);
    end
    load(calibration_file, 'Calib', 'SS0');
    fprintf('\n=== Loaded explicit E-TANK calibration: %s ===\n', calibration_file);
else
    fprintf('\n=== Using embedded final thesis calibration ===\n');

    % Final calibrated structural parameters and policy normalizations.
    Calib = struct();
    Calib.beta = 0.997;
    Calib.sigma = 1.5;
    Calib.psi = 2;
    Calib.chi = 36.984023760240923;
    Calib.lambda_htM = 0.20000000000000001;
    Calib.omega = 0.012;
    Calib.epsilon_E = 0.29999999999999999;
    Calib.ebar = 0.048810491265711062;
    Calib.alpha_y = 0.29999999999999999;
    Calib.beta_y = 0.63;
    Calib.alpha_c = 0.20000000000000001;
    Calib.sigma_L = 1.3999999999999999;
    Calib.omega_L = 0.81912080591756276;
    Calib.delta = 0.025000000000000001;
    Calib.phi_k = 4;
    Calib.ass = 1.7559813768127772;
    Calib.A = 1.7559813768127772;
    Calib.alpha_E = 0.25;
    Calib.rho_E = 1.8;
    Calib.phi_mix = 75.410145028683701;
    Calib.gz = 1.0029999999999999;
    Calib.delta_c = 0.0050000000000000001;
    Calib.eta_c = 7.7217487048955977;
    Calib.phi_d = 1;
    Calib.delta_clim = 0.0035000000000000001;
    Calib.d0 = 0;
    Calib.d1 = 0;
    Calib.d2 = 1.6023073225444642e-08;
    Calib.eta_X = 0.59999999999999998;
    Calib.D_ex_ss = 0.13696323739659313;
    Calib.gz_ss = 1.0029999999999999;
    Calib.tau_ss = 0.050000000000000003;
    Calib.P_oil_ss = 1;
    Calib.A_c_ss = 1.0000000000000024;
    Calib.alpha_KL = 0.32258064516129031;
    Calib.m_Y = 0.99562965575655393;
    Calib.sigma_Y = 0.5;
    Calib.theta0 = 0.20000000000000001;
    Calib.phi0 = 0.80000000000000004;
    Calib.tau0 = 0.050000000000000003;

    % Complete common initial steady state associated with Calib.
    SS0 = struct();
    SS0.P_Z = 1.0077629025699417;
    SS0.Z_H = 0.41608728229586334;
    SS0.Z_R = 0.65878666076902426;
    SS0.C_H_n = 0.41204903040804042;
    SS0.C_R_n = 0.65239293861091363;
    SS0.C_H_e = 0.053075492317526181;
    SS0.C_R_e = 0.055563223013719097;
    SS0.L_H_u = 0.36352999479065101;
    SS0.L_R_s = 0.32578416796900161;
    SS0.Lambda_R = 1.8557712293223738;
    SS0.I = 0.25871260563308529;
    SS0.K = 9.2337276293386221;
    SS0.q_K = 1.0122204127236445;
    SS0.E = 0.096141334195180159;
    SS0.E_c = 0.019225351297655697;
    SS0.E_d = 0.13696323739659483;
    SS0.E_y = 0.041075657320699649;
    SS0.m = 0.14036869793013287;
    SS0.P_e = 1.7041723630488299;
    SS0.P_c = 1.0418538057136237;
    SS0.P_d = 1.0499999999999923;
    SS0.er = 0.99999999999999234;
    SS0.Y_g = 1.1552788451976113;
    SS0.Y_net = 1.0000000000000178;
    SS0.r_k = 0.0329234317154002;
    SS0.L_a = 0.20813103977352468;
    SS0.L_s_y = 0.25997981265556414;
    SS0.L_u_y = 0.06058491766584876;
    SS0.L_u_c = 0.012121081292281432;
    SS0.K_c = 0.12167629175955401;
    SS0.L_R_c = 0.00064752171963715875;
    SS0.w_s = 2.1151914019709288;
    SS0.w_u = 1.3219946262318247;
    SS0.W_a = 3.0269391854551766;
    SS0.V_KL = 0.7043390879766579;
    SS0.A_c = 1.0000000000000024;
    SS0.gamma_agg = 1.0029999999999999;
    SS0.e_em = 0.13696323739659483;
    SS0.X = 2896.2752106847151;
    SS0.Dam = 0.13440810921369639;
    SS0.R_tau = 0.0068481618698297414;
    SS0.T_H = 0.021914117983455175;
    SS0.T_R = 0.001369632373965948;
    SS0.C_n = 0.60432415697033903;
    SS0.EX = 0.13696323739659377;
    SS0.W_H = -822.05851411057142;
    SS0.W_R = -643.4264481935943;
    SS0.tau = 0.050000000000000003;
    SS0.P_oil = 1;
    SS0.e_ROW = 10;
    SS0.theta_policy = 0.20000000000000001;
    SS0.phi_policy = 0.80000000000000004;
    SS0.burden_H = 0.1799999999999968;
    SS0.burden_R = 0.12674549469868859;
    SS0.burden_ratio = 1.4201688227886118;
end

% Make every calibrated scalar available to the MATLAB and Dynare blocks.
calibration_fields = fieldnames(Calib);
for calibration_i = 1:numel(calibration_fields)
    calibration_name = calibration_fields{calibration_i};
    eval([calibration_name ' = Calib.(calibration_name);']);
end
steady_state_fields = fieldnames(SS0);
for calibration_i = 1:numel(steady_state_fields)
    calibration_name = steady_state_fields{calibration_i};
    if isnumeric(SS0.(calibration_name)) && isscalar(SS0.(calibration_name))
        eval([calibration_name '_ss = SS0.(calibration_name);']);
    end
end
theta_benchmark = Calib.theta0;
phi_benchmark = Calib.phi0;
K_total_ss = SS0.K;

% Loading or embedding MATLAB variables is not sufficient for Dynare:
% synchronize every structural parameter explicitly with M_.params.
parameter_names = cellstr(M_.param_names);
for calibration_i = 1:numel(parameter_names)
    calibration_name = strtrim(parameter_names{calibration_i});
    if isfield(Calib, calibration_name)
        set_param_value(calibration_name, Calib.(calibration_name));
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

% [9] Ricardian budget constraint
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

% [17] Clean-energy price: Cobb-Douglas unit cost
% P_c = (1/A_c)*(r_k/alpha_c)^alpha_c * (w_u/(1-alpha_c))^(1-alpha_c)
P_c = (1/A_c) * (r_k/alpha_c)^alpha_c * (w_u/(1-alpha_c))^(1-alpha_c);

% [18] Retailer FOC for clean energy (with mix friction)
% P_e*alpha_E*(E/E_c)^(1/rho_E) = P_c + phi_mix*(m - m(-1))*P_e*E/E_d
P_e*alpha_E*(E/E_c)^(1/rho_E) = P_c + phi_mix*(m - m(-1))*P_e*E/E_d;

% [19] Retailer FOC for dirty energy (with mix friction)
% P_e*(1-alpha_E)*(E/E_d)^(1/rho_E) = P_d - phi_mix*(m - m(-1))*E_c*P_e*E/E_d^2
P_e*(1-alpha_E)*(E/E_d)^(1/rho_E) = P_d - phi_mix*(m - m(-1))*E_c*P_e*E/E_d^2;

% =========================================================================
% BLOCK 3: CAPITAL-LABOR-ENERGY FINAL PRODUCTION
% =========================================================================

% [20] Capital-labor composite (inner Cobb-Douglas, constant returns)
V_KL = (K(-1)-K_c)^alpha_KL * L_a^(1-alpha_KL);

% [20b] nested CES final production : composite KL + energy
Y_g = ass * ( m_Y*V_KL^(1-1/sigma_Y) + (1-m_Y)*E_y^(1-1/sigma_Y) )^(1/(1-1/sigma_Y));

% [21] Output net of damages
Y_net = (1 - Dam) * Y_g;

% [22] Capital FOC (outer CES, via the V_KL composite)
% The factor ass^(1-1/sigma_Y) is required by the standard CES derivative
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

% [28] Unskilled-labor factor demand in clean energy
% L_u_c = (1-alpha_c)*P_c*E_c/w_u  (labor share of cost, zero profit under CRS)
L_u_c = (1-alpha_c)*P_c*E_c/w_u;

% [46] Capital factor demand in clean energy
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
% BLOCK 5: REDUCED-FORM QUADRATIC CLIMATE DAMAGES
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
% theta_policy and phi_policy are deterministic exogenous policy variables.
% Period t=0 is the common pre-policy baseline; the alternative recycling
% rule applies from t=1 onward. Structural parameters remain fixed.

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
% Numerical solution horizon: 1,600 quarters (400 years) by default.
% This long horizon limits terminal-boundary effects; it is not a forecast.
% Additive carbon-price wedge: 30-year ramp from tau_init to the terminal value

@#ifndef T_SIM_QUARTERS
@#define T_SIM_QUARTERS = 1600
@#endif
@#define T_sim = T_SIM_QUARTERS
@#ifndef TAU_FINAL
@#define TAU_FINAL = 2.0
@#endif
@#ifndef KAPPA_TAX
@#define KAPPA_TAX = 0.5
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
K        = K_total_ss;   % K_total = K_y + K_c
K_c      = K_c_ss;       % initial steady-state clean-energy capital
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

% Terminal conditions. TAU_FINAL is used consistently by endval and by the
% full exogenous trajectory below.
endval;
tau      = @{TAU_FINAL};
P_oil    = P_oil_ss;
e_ROW    = e_ROW_ss;
theta_policy = @{THETA_POLICY};
phi_policy   = @{PHI_POLICY};
end;

% Policy path: tau_init to tau_final over 120 quarters (power-function ramp)
% The exact trajectory is injected further below, directly into this .mod,
% into oo_.exo_simul after perfect_foresight_setup.
shocks;
var tau;
periods 1;
values 0.05;
end;
% The block below, after perfect_foresight_setup, overrides
% oo_.exo_simul with the exact trajectory and configured curvature.

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

%% Inject the configured power-function tau trajectory into oo_.exo_simul
% Override the shocks-block placeholder with the complete transition path
tau_idx_v   = find(strcmp(M_.exo_names, 'tau'));
poil_idx_v  = find(strcmp(M_.exo_names, 'P_oil'));
tau_init_v  = 0.05;
tau_final_v = @{TAU_FINAL};
T_carbon_v  = 120;
kappa_v     = @{KAPPA_TAX};
T_sim_v     = @{T_sim};   % simulation horizon in quarters
% oo_.exo_simul is (T_sim+2) x n_exo:
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

    % --- Household energy burdens ---
    burden_H_sim = P_e_sim.*C_H_e_sim ./ (C_H_n_sim + P_e_sim.*C_H_e_sim);
    burden_R_sim = P_e_sim.*C_R_e_sim ./ (C_R_n_sim + P_e_sim.*C_R_e_sim);
    regress_sim  = burden_H_sim ./ burden_R_sim;
    prem_sal_sim = w_s_sim ./ w_u_sim;

    % Normalize against the common pre-announcement baseline,
    % not the terminal policy-specific steady state or traj(t=1).
    Y_ss   = SS0.Y_net;
    Pe_ss  = SS0.P_e;
    Ac_ss  = SS0.A_c;
    WH_ss  = SS0.W_H;
    WR_ss  = SS0.W_R;

    emissions_reduction_30 = 100*(1-e_em_sim(120)/SS0.e_em);
    Y_net_deviation_30 = 100*(Y_net_sim(120)/Y_ss-1);
    A_c_30 = A_c_sim(120);
    clean_share_30 = 100*E_c_sim(120)/(E_c_sim(120)+E_d_sim(120));

    fprintf('\n=== RESULTS AT 30 YEARS (quarter 120) ===\n');
    fprintf('Domestic emissions reduction : %.15g%%\n', emissions_reduction_30);
    fprintf('Y_net deviation              : %+.15g%% vs SS\n', Y_net_deviation_30);
    fprintf('A_c level                    : %.15g\n', A_c_30);
    fprintf('Clean-energy share           : %.15g%%\n', clean_share_30);
    fprintf('P_e       : %+.2f%% vs SS\n', 100*(P_e_sim(120)/Pe_ss-1));
    fprintf('burden_H  : %.4f (SS: %.4f)\n', burden_H_sim(120), P_e_ss*C_H_e_ss/(C_H_n_ss+P_e_ss*C_H_e_ss));
    fprintf('burden ratio: %.4f (SS: %.4f)\n', regress_sim(120), SS0.burden_ratio);

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
    % Convert type-specific S-CEVs to a common normalized money metric before
    % population aggregation, exactly as in the final production postprocessor.
    money_scale_H = SS0.P_Z*SS0.Z_H/(SS0.C_H_n + SS0.P_e*SS0.C_H_e);
    money_scale_R = SS0.P_Z*SS0.Z_R/(SS0.C_R_n + SS0.P_e*SS0.C_R_e);
    money_metric_H = SCEV_H*money_scale_H;
    money_metric_R = SCEV_R*money_scale_R;
    money_metric_aggregate = lambda_htM*money_metric_H + (1-lambda_htM)*money_metric_R;
    initial_transfer_ratio = SS0.T_H/SS0.T_R;

    fprintf('S-CEV HtM                  : %+.15g%%\n', SCEV_H);
    fprintf('S-CEV Ricardian            : %+.15g%%\n', SCEV_R);
    fprintf('Aggregate money metric     : %+.15g\n', money_metric_aggregate);
    fprintf('Initial burdens H/R/ratio  : %.15g / %.15g / %.15g\n', ...
        SS0.burden_H, SS0.burden_R, SS0.burden_ratio);
    fprintf('Initial transfer ratio H/R : %.15g\n', initial_transfer_ratio);

else
    warning('Simulation did not converge. Check initval/endval and the residuals.');
end
