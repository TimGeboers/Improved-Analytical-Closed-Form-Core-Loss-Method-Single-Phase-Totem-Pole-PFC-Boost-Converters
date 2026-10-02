%% Core Loss Calculation for Coupled Inductors using Analytical iGSE
% ===========================================================================================================
% Paper Reference: "Improved Core Loss Calculation Method for Single-Phase Totem-Pole PFC Boost Converters"
% Paper Authors: Tim Geboers, Wout Vanderwegen, Wilmar Martinez, Camilo Suarez
% ===========================================================================================================
% Code Authors: Tim Geboers, Camilo Suarez
% ===========================================================================================================
% (c) 2026, Tim Geboers
% ===========================================================================================================

clear; clc; close all;
tic; % Start execution timer

% =======================================
%%  1. CONVERTER & ELECTRICAL PARAMETERS  
% =======================================
Vin_rms = 120;                  % V (Grid input RMS)
Vi_pk   = sqrt(2) * Vin_rms;    % V (Grid input peak)
Vo      = 400;                  % V (DC Output voltage)
Po      = 3300;                 % W (Total output power)

Iin_rms = Po / Vin_rms;         % A (Input RMS current)
Iin_pk  = sqrt(2) * Iin_rms;    % A (Input peak current)
Iph_pk  = Iin_pk / 2;           % A (Per-phase peak current)

fp      = 50;                   % Hz (Line frequency)
Tp      = 1 / fp;               % s  (Line period)
omega_p = 2 * pi / Tp;          % rad/s (Line angular frequency)

fsw     = 120e3;                % Hz (Switching frequency)
Ts      = 1 / fsw;              % s  (Switching period)

% =============================================================
%%  2. MAGNETIC CORE GEOMETRY & MATERIAL PROPERTIES (SI Units)
% =============================================================
% Geometry dimensions (3-leg C-core structure)
w_window = 10e-3;               % m (Window width)
w_core   = 17e-3;               % m (Outer leg width)
h_window = 11e-3;               % m (Window height)
d_core   = 29.14e-3;            % m (Core depth)
h_midleg = 13.36e-3;            % m (Center leg height)
l_pole   = 4.88e-3;             % m (Pole length)
lg_mag   = 0.1e-3;              % m (Air gap length)

% Derived core path lengths and cross-sectional areas
l_c_core  = (w_window + (2*w_core) + (2*h_window) - (2*lg_mag) + h_midleg); % m (Core path length)
Ae_c_core = w_core * d_core;    % m^2 (Outer leg cross-sectional area, Ae,out)
Ae_T_core = h_midleg * d_core;  % m^2 (Center leg cross-sectional area, Ae,ctr)

% Derived core volumes
Ve_c_core = l_c_core * Ae_c_core;   % m^3 (Outer leg volume, Ve,out)
V_pole    = Ae_T_core * l_pole;     % m^3 (Single pole volume)
Ve_T_core = 2 * V_pole;             % m^3 (Center leg total volume, Ve,ctr)

% Winding & Inductance properties (Section II-B)
N   = 6;                        % Number of turns per phase
L11 = 66.1e-6;                  % H (Self-inductance Phase 1)
L12 = 37.6e-6;                  
M   = L12;                      % H (Mutual inductance)
Lleak = L11 - M;                % H (Equivalent leakage inductance, Eq. 5)

% iGSE / Steinmetz Material Coefficients (Section III-A)
alpha = 1.6931;
beta  = 2.8325;
kc    = 0.1172;

% Calculate iGSE constant ki (Section III-B)
theta_int = linspace(0, 2*pi, 1000);
integral_iGSE = trapz(theta_int, abs(cos(theta_int)).^alpha);
ki_iGSE = kc / ( 2^(beta - alpha) * (2*pi)^(alpha - 1) * integral_iGSE );

% ===========================================================
%%  3. TIME-DOMAIN WAVEFORMS & FLUX DERIVATIONS (Section II)
% ===========================================================
% Phase angle t0 calculation for CCM integration limits [s]
t0 = Tp/pi * atan( ( sqrt(4*pi^2*Lleak^2*(Iph_pk)^2 + Tp^2*Vi_pk^2) - Tp*Vi_pk ) ...
                   / ( 2*pi*Lleak*Iph_pk ) );
theta0 = omega_p * t0; % [rad]

%% --- 3.1 Rising Portion (t = to .. Tp/4) ---
t_r = (t0 + Ts) : Ts : (Tp/4); % Rising time period [s]

% Solve Eq. 3 - Grid input voltage [V]
vin_r = Vi_pk * sin(omega_p * t_r);   

% Solve Eq. 4 - Duty cycle vectorization [-]
d_r = ( Vo - ( vin_r - omega_p*Lleak*Iph_pk.*cos(omega_p*t_r) ) ) / Vo;
d_r = min(max(d_r, 0), 1); % Enforce CCM bounds [0, 1]

% Solve Eq. 1 - Average phase current [A]
iavg_r = Po * Vo * (1 - d_r) / (2 * Vin_rms^2);

% Solve Eq. 2 - Flux peak-to-peak ripple amplitude [Wb]
delta_phi_r = (Vo - vin_r) .* (1 - d_r) * Ts / N;

% Solve Eq. 6 & 7 - Outer leg flux upper and lower peak envelope + Outer leg flux [Wb]
phi_high_r = Lleak .* iavg_r / N + 0.5 * delta_phi_r;
phi_low_r = Lleak .* iavg_r / N - 0.5 * delta_phi_r;
phi_tot_r_outer = reshape([phi_high_r; phi_low_r], 1, []);
t_phi_r = reshape([t_r; t_r + (1 - d_r)*Ts], 1, []); % [s]

% Solve Eq. 8 - Center leg flux [Wb]
phi_high_r_ctr = (2 * Lleak .* iavg_r / N) + 0.5 * abs(delta_phi_r .* (1 - 2*d_r));
phi_low_r_ctr = (2 * Lleak .* iavg_r / N) - 0.5 * abs(delta_phi_r .* (1 - 2*d_r));
phi_tot_r_ctr = reshape([phi_high_r_ctr; phi_low_r_ctr], 1, []);

% Solve Eq. 13 - Major loop peak-to-peak flux density excursion [T]
Delta_BM_r_outer = 2 * max(phi_tot_r_outer) / Ae_c_core; % Outer leg
Delta_BM_r_ctr = 2 * max(phi_tot_r_ctr) / Ae_T_core; % Center leg

%% --- 3.2 Falling Portion (t = Tp/4 .. Tp/2 - to) ---
t_f = (Tp/4) : Ts : (Tp/2 - t0); % Falling time period [s]

% Solve Eq. 3 - Grid input voltage [V]
vin_f = Vi_pk * sin(omega_p * t_f);

% Solve Eq. 4 - Duty cycle vectorization [-]
d_f = ( Vo - ( vin_f - omega_p*Lleak*Iph_pk.*cos(omega_p*t_f) ) ) / Vo;
d_f = min(max(d_f, 0), 1); % Enforce CCM bounds [0, 1]

% Solve Eq. 1 - Average phase current [A]
iavg_f = Po * Vo * (1 - d_f) / (2 * Vin_rms^2);

% Solve Eq. 2 - Flux peak-to-peak ripple amplitude [Wb]
delta_phi_f = (Vo - vin_f) .* (1 - d_f) * Ts / N;   

% Solve Eq. 6 & 7 - Outer leg flux upper and lower peak envelope + Outer leg flux [Wb]
phi_high_f = Lleak .* iavg_f / N + 0.5 * delta_phi_f;
phi_low_f = Lleak .* iavg_f / N - 0.5 * delta_phi_f;
phi_tot_f_outer = reshape([phi_high_f; phi_low_f], 1, []);

% Solve Eq. 13 - Major loop peak-to-peak flux density excursion [T]
Delta_BM_f_outer = Delta_BM_r_outer; % Outer leg (Same Delta_B for rising and falling portion)
Delta_BM_f_ctr = Delta_BM_r_ctr; % Center leg (Same Delta_B for rising and falling portion)

% ===================================================
%%  4. ANALYTICAL CORE LOSS EVALUATION (Section III)
% ===================================================

% --------------------------------------------------------------------------------
%% 4.1 MAJOR LOOP LOSSES - ANALYTICAL INTEGRATION (Section III-B, Eq. 9 & 10 & 27)
% --------------------------------------------------------------------------------
%% --- Outer Leg Major Losses ---
% Solve Eq. 9 & 10 & 27 - Volumetric iGSE core loss (integral) [W/m^3]
% --- Rising Portion ---
pref = (ki_iGSE * (Delta_BM_r_outer)^(beta - alpha)) / (Tp * (N * Ae_c_core)^alpha);
integrand_vec = (abs(vin_r).^alpha) .* (d_r - ((Vo - vin_r) .* (1 - d_r)) ./ (vin_r + eps));
PM_r_Q1 = pref * trapz(t_r, integrand_vec);
PM_r_outer_integral = 2 * PM_r_Q1; % accounts for two quarter periods (Q1 + Q3) for major rising portion

% --- Falling Portion ---
pref_falling = (ki_iGSE * (Delta_BM_f_outer)^(beta - alpha)) / (Tp * (N * Ae_c_core)^alpha);
weight_f = (1 - d_f - (vin_f .* d_f) ./ (Vo - vin_f + eps));
integrand_vec_falling = (abs(Vo - vin_f).^alpha) .* weight_f;
PM_f_Q2 = pref_falling * trapz(t_f, integrand_vec_falling);
PM_f_outer_integral = 2 * PM_f_Q2; % accounts for two quarter periods (Q2 + Q4) for major falling portion

% --- Rising + Falling Portion ---
PM_outer_vol_integral  = PM_r_outer_integral + PM_f_outer_integral; 

%% --- Center Leg Major Losses ---
% Solve Eq. 9 & 10 & 27 - Volumetric iGSE core loss (integral) [W/m^3]
% --- Rising Portion ---
pref_central = (ki_iGSE * (Delta_BM_r_ctr)^(beta - alpha)) / (Tp * (N * Ae_T_core)^alpha);
integrand_vec_central = (abs(vin_r).^alpha) .* (d_r - ((Vo - vin_r) .* (1 - d_r)) ./ (vin_r + eps));
PM_r_Q1_central = pref_central * trapz(t_r, integrand_vec_central);
PM_r_ctr_integral = 2 * PM_r_Q1_central; % accounts for two quarter periods (Q1 + Q3) for major rising portion

% --- Falling Portion ---
pref_falling_central = (ki_iGSE * (Delta_BM_f_ctr)^(beta - alpha)) / (Tp * (N * Ae_T_core)^alpha);
weight_f_central = (1 - d_f - (vin_f .* d_f) ./ (Vo - vin_f + eps));
integrand_vec_falling_central = (abs(Vo - vin_f).^alpha) .* weight_f_central;
PM_f_Q2_central = pref_falling_central * trapz(t_f, integrand_vec_falling_central);
PM_f_ctr_integral = 2 * PM_f_Q2_central; % accounts for two quarter periods (Q2 + Q4) for major falling portion

% --- Rising + Falling Portion ---
PM_ctr_vol_integral  = PM_r_ctr_integral + PM_f_ctr_integral; 

% ---------------------------------------------------------------------------------
%% 4.2 MINOR LOOP LOSSES - ANALYTICAL INTEGRATION (Section III-C, Eq. 16 & 17 & 28)
% ---------------------------------------------------------------------------------
%% --- Outer Leg Minor Losses ---
% Solve Eq. 16 & 17 & 28 - Volumetric iGSE core loss (integral) [W/m^3]
% --- Rising Portion ---
dB_ripple = (Ts * (Vo - vin_r) .* (1 - d_r)) ./ (N * Ae_c_core);
pref_minor = (ki_iGSE) / (Tp * (N * Ae_c_core)^alpha);
integrand = ki_iGSE * (abs(dB_ripple)).^(beta - alpha) .* ...
            ( (abs(Vo - vin_r) / (N*Ae_c_core)).^alpha .* (1 - d_r) + ...
              (abs(vin_r)      / (N*Ae_c_core)).^alpha .* (d_r) );
PI_outer_leg_integral = (1/Tp) * trapz(t_r, integrand);
Pm_r_outer_integral = 2 * PI_outer_leg_integral; % accounts for two quarter periods (Q1 + Q3) for minor rising portion

% --- Falling Portion ---
dB_ripple_fall = (Ts * vin_f .* d_f) ./ (N * Ae_c_core);
integrand_fall = ki_iGSE * (abs(dB_ripple_fall)).^(beta - alpha) .* ...
                 ( (abs(Vo - vin_f) / (N*Ae_c_core)).^alpha .* (1 - d_f) + ...
                   (abs(vin_f)      / (N*Ae_c_core)).^alpha .* (d_f) );
PI_falling_outer_leg_integral = (1/Tp) * trapz(t_f, integrand_fall);
Pm_f_outer_integral = 2 * PI_falling_outer_leg_integral; % accounts for two quarter periods (Q2 + Q4) for minor falling portion

% --- Rising + Falling Portion ---
Pm_outer_vol_integral  = Pm_r_outer_integral + Pm_f_outer_integral; 

%% --- Center Leg Minor Losses ---
% Solve Eq. 16 & 17 & 28 - Volumetric iGSE core loss (integral) [W/m^3]
% --- Rising Portion ---
dB_ripple_central = (Ts * (Vo - vin_r) .* (1 - d_r)) ./ (N * Ae_T_core);
pref_minor_central = (ki_iGSE) / (Tp * (N * Ae_T_core)^alpha);
integrand_central = ki_iGSE * (abs(dB_ripple_central)).^(beta - alpha) .* ...
            ( (abs(Vo - vin_r) / (N*Ae_T_core)).^alpha .* (1 - d_r) + ...
              (abs(vin_r)      / (N*Ae_T_core)).^alpha .* (d_r) );
PI_central_leg_integral = (1/Tp) * trapz(t_r, integrand_central);
Pm_r_ctr_integral = 2 * PI_central_leg_integral; % accounts for two quarter periods (Q1 + Q3) for minor rising portion

% --- Falling Portion ---
dB_ripple_fall_central = (Ts * vin_f .* d_f) ./ (N * Ae_T_core);
integrand_fall_central = ki_iGSE * (abs(dB_ripple_fall_central)).^(beta - alpha) .* ...
                 ( (abs(Vo - vin_f) / (N*Ae_T_core)).^alpha .* (1 - d_f) + ...
                   (abs(vin_f)      / (N*Ae_T_core)).^alpha .* (d_f) );
PI_falling_central_leg_integral = (1/Tp) * trapz(t_f, integrand_fall_central);
Pm_f_ctr_integral = 2 * PI_falling_central_leg_integral; % accounts for two quarter periods (Q2 + Q4) for minor falling portion

% --- Rising + Falling Portion ---
Pm_ctr_vol_integral  = Pm_r_ctr_integral + Pm_f_ctr_integral; 

% ----------------------------------------------------------------------
%% 4.3 MAJOR LOOP LOSSES - CLOSED FORM (Section III-B, Eq. 14 & 15 & 27)
% ----------------------------------------------------------------------
%% --- Outer Leg Major Losses ---
% Solve Eq. 14 & 15 & 27 - Volumetric iGSE core loss (closed form) [W/m^3]
% --- Rising Portion ---
pref_M_r_outer = (2 * ki_iGSE * (Delta_BM_r_outer)^(beta - alpha)) / (Tp * (N * Ae_c_core)^alpha);
PM_r_outer_CF  = pref_M_r_outer * ( (Lleak * Iph_pk * Vi_pk^(alpha - 1) / alpha) * (1 - (sin(theta0))^alpha) );

% --- Falling Portion ---
pref_M_f_outer = (2 * ki_iGSE * (Delta_BM_f_outer)^(beta - alpha)) / (Tp * (N * Ae_c_core)^alpha);
term_upper_out = (Vo - Vi_pk * sin(theta0))^alpha;
term_lower_out = (Vo - Vi_pk)^alpha;
PM_f_outer_CF  = pref_M_f_outer * ((Lleak * Iph_pk) / (Vi_pk * alpha)) * (term_upper_out - term_lower_out);

% --- Rising + Falling Portion ---
PM_outer_vol_CF = PM_r_outer_CF + PM_f_outer_CF;

%% --- Center Leg Major Losses ---
% Solve Eq. 14 & 15 & 27 - Volumetric iGSE core loss (closed form) [W/m^3]
% --- Rising Portion ---
pref_M_r_ctr = (2 * ki_iGSE * (Delta_BM_r_ctr)^(beta - alpha)) / (Tp * (N * Ae_T_core)^alpha);
PM_r_ctr_CF  = pref_M_r_ctr * ( (Lleak * Iph_pk * Vi_pk^(alpha - 1) / alpha) * (1 - (sin(theta0))^alpha) );

% --- Falling Portion ---
pref_M_f_ctr = (2 * ki_iGSE * (Delta_BM_f_ctr)^(beta - alpha)) / (Tp * (N * Ae_T_core)^alpha);
PM_f_ctr_CF  = pref_M_f_ctr * ((Lleak * Iph_pk) / (Vi_pk * alpha)) * (term_upper_out - term_lower_out);

% --- Rising + Falling Portion ---
PM_ctr_vol_CF   = PM_r_ctr_CF + PM_f_ctr_CF;

% ---------------------------------------------------------------------------------------
%% 4.4 MINOR LOOP LOSSES - CLOSED FORM TAYLOR EXPANSION (Section III-C, Eq. 25 & 26 & 28)
% ---------------------------------------------------------------------------------------
max_order = 12; % Adapt this value (maximum order K) depending on the desired accuracy 
num_terms = max_order + 1;
syms theta

M_mod = Vi_pk / Vo;
gamma = (2 * pi / Tp) * Lleak * Iph_pk / Vo;

%% --- Outer Leg Minor Losses ---
% --- Symbolic Integrands in ANGULAR DOMAIN for Outer Leg ---
dBm_r_out = (Ts * Vo / (N * Ae_c_core)) * (1 - M_mod*sin(theta)) * (M_mod*sin(theta) - gamma*cos(theta));
dBm_f_out = (Ts * Vo / (N * Ae_c_core)) * (M_mod*sin(theta)) * (1 - M_mod*sin(theta) + gamma*cos(theta));
um_x_out = (Vo * (1 - M_mod*sin(theta)))^alpha * (M_mod*sin(theta) - gamma*cos(theta)) + ...
             (Vo * M_mod*sin(theta))^alpha * (1 - M_mod*sin(theta) + gamma*cos(theta));

% Solve Eq. 22 & 23 - Generalized algebraic polynomials [-]
integrand_r_out = (dBm_r_out)^(beta - alpha) * um_x_out;
integrand_f_out = (dBm_f_out)^(beta - alpha) * um_x_out;

theta_r = pi / 4;   % Rising Phase Interval (Section III-C)
theta_f = 3 * pi / 4; % Falling Phase Interval (Section III-C)

% Solve Eq. 24 - Scalar coefficients [-]
A_out = zeros(1, num_terms);
B_out = zeros(1, num_terms);
for n = 0:max_order
    A_out(n+1) = double(subs(diff(integrand_r_out, theta, n), theta, theta_r) / factorial(n)); 
    B_out(n+1) = double(subs(diff(integrand_f_out, theta, n), theta, theta_f) / factorial(n)); 
end

% Solve Eq. 25 & 26 & 28 - Volumetric iGSE core loss (closed form) [W/m^3]
front_mult_out = ki_iGSE / (2 * pi * (N * Ae_c_core)^alpha);

% --- Rising Portion ---
Pm_r_int_out = sum( (A_out ./ (1:num_terms)) .* ((pi/2 - theta_r).^(1:num_terms) - (theta0 - theta_r).^(1:num_terms)) );
Pm_r_outer_CF = 2 * front_mult_out * Pm_r_int_out; % accounts for two quarter periods (Q1 + Q3) for minor rising portion

% --- Falling Portion ---
Pm_f_int_out = sum( (B_out ./ (1:num_terms)) .* (((pi - theta0) - theta_f).^(1:num_terms) - (pi/2 - theta_f).^(1:num_terms)) );
Pm_f_outer_CF = 2 * front_mult_out * Pm_f_int_out; % accounts for two quarter periods (Q2 + Q4) for minor falling portion

% --- Rising + Falling Portion ---
Pm_outer_vol_CF  = Pm_r_outer_CF + Pm_f_outer_CF; 

%% --- Center Leg Minor Losses ---
% --- Symbolic Integrands in ANGULAR DOMAIN for Center Leg ---
dBm_r_ctr = (Ts * Vo / (N * Ae_T_core)) * (1 - M_mod*sin(theta)) * (M_mod*sin(theta) - gamma*cos(theta));
dBm_f_ctr = (Ts * Vo / (N * Ae_T_core)) * (M_mod*sin(theta)) * (1 - M_mod*sin(theta) + gamma*cos(theta));
um_x_ctr = (Vo * (1 - M_mod*sin(theta)))^alpha * (M_mod*sin(theta) - gamma*cos(theta)) + ...
             (Vo * M_mod*sin(theta))^alpha * (1 - M_mod*sin(theta) + gamma*cos(theta));

% Solve Eq. 22 & 23 - Generalized algebraic polynomials [-]
integrand_r_ctr = (dBm_r_ctr)^(beta - alpha) * um_x_ctr;
integrand_f_ctr = (dBm_f_ctr)^(beta - alpha) * um_x_ctr;

% Solve Eq. 24 - Scalar coefficients [-]
A_ctr = zeros(1, num_terms);
B_ctr = zeros(1, num_terms);
for n = 0:max_order
    A_ctr(n+1) = double(subs(diff(integrand_r_ctr, theta, n), theta, theta_r) / factorial(n)); 
    B_ctr(n+1) = double(subs(diff(integrand_f_ctr, theta, n), theta, theta_f) / factorial(n)); 
end

% Solve Eq. 25 & 26 & 28 - Volumetric iGSE core loss (closed form) [W/m^3]
front_mult_ctr = ki_iGSE / (2 * pi * (N * Ae_T_core)^alpha);

% --- Rising Portion ---
Pm_r_int_ctr = sum( (A_ctr ./ (1:num_terms)) .* ((pi/2 - theta_r).^(1:num_terms) - (theta0 - theta_r).^(1:num_terms)) );
Pm_r_ctr_CF = 2 * front_mult_ctr * Pm_r_int_ctr; % accounts for two quarter periods (Q1 + Q3) for minor rising portion

% --- Falling Portion ---
Pm_f_int_ctr = sum( (B_ctr ./ (1:num_terms)) .* (((pi - theta0) - theta_f).^(1:num_terms) - (pi/2 - theta_f).^(1:num_terms)) );
Pm_f_ctr_CF = 2 * front_mult_ctr * Pm_f_int_ctr; % accounts for two quarter periods (Q2 + Q4) for minor falling portion

% --- Rising + Falling Portion ---
Pm_ctr_vol_CF  = Pm_r_ctr_CF + Pm_f_ctr_CF;

% ========================================================
%%  5. TOTAL LOSS COMPUTATION & REPORTING (Section III-D)
% ========================================================
% Closed Form Power Loss per Core Segment [W]
P_outer_legs_CF = 2 * (PM_outer_vol_CF + Pm_outer_vol_CF) * Ve_c_core; % 2 Outer Legs
P_center_leg_CF = (PM_ctr_vol_CF + Pm_ctr_vol_CF) * Ve_T_core; % 1 Center Leg

% Total Closed Form Power Loss [W]
PM_total_CF = (2 * PM_outer_vol_CF * Ve_c_core) + (PM_ctr_vol_CF * Ve_T_core); % W (Total Major Loop Core Loss)
Pm_total_CF = (2 * Pm_outer_vol_CF * Ve_c_core) + (Pm_ctr_vol_CF * Ve_T_core); % W (Total Minor Loop Core Loss)
P_total_core_CF = P_outer_legs_CF + P_center_leg_CF; % W (Total Core Loss)

% Integral Power Loss per Core Segment [W]
P_outer_legs_integral = 2 * (PM_outer_vol_integral + Pm_outer_vol_integral) * Ve_c_core; % 2 Outer Legs
P_center_leg_integral = (PM_ctr_vol_integral + Pm_ctr_vol_integral) * Ve_T_core; % 1 Center Leg

% Total Integral Power Loss [W]
PM_total_integral = (2 * PM_outer_vol_integral * Ve_c_core) + (PM_ctr_vol_integral * Ve_T_core); % W (Total Major Loop Core Loss)
Pm_total_integral = (2 * Pm_outer_vol_integral * Ve_c_core) + (Pm_ctr_vol_integral * Ve_T_core); % W (Total Minor Loop Core Loss)
P_total_core_integral = P_outer_legs_integral + P_center_leg_integral; % W (Total Core Loss)

% Display Results in Command Window
fprintf('\n=============================================================\n');
fprintf('   ANALYTICAL CLOSED FORM iGSE CORE LOSS CALCULATION REPORT\n');
fprintf('=============================================================\n');
fprintf('Major Loop Loss (2 Outer Legs):       %.4f W\n', 2 * PM_outer_vol_CF * Ve_c_core);
fprintf('Minor Loop Loss (2 Outer Legs):       %.4f W\n', 2 * Pm_outer_vol_CF * Ve_c_core);
fprintf('Total Loss (2 Outer Legs):          %.4f W\n', P_outer_legs_CF);
fprintf('-------------------------------------------------------------\n');
fprintf('Major Loop Loss (Center Leg):       %.4f W\n', PM_ctr_vol_CF * Ve_T_core);
fprintf('Minor Loop Loss (Center Leg):       %.4f W\n', Pm_ctr_vol_CF * Ve_T_core);
fprintf('Total Loss (Center Leg):            %.4f W\n', P_center_leg_CF);
fprintf('=============================================================\n');
fprintf('Major Loop Loss (Total):       %.4f W\n', PM_total_CF);
fprintf('Minor Loop Loss (Total):       %.4f W\n', Pm_total_CF);
fprintf('TOTAL COUPLED INDUCTOR CORE LOSS (closed form):   %.4f W\n', P_total_core_CF);
fprintf('=============================================================\n');

% Display Results in Command Window
fprintf('\n=============================================================\n');
fprintf('   ANALYTICAL INTEGRAL iGSE CORE LOSS CALCULATION REPORT\n');
fprintf('=============================================================\n');
fprintf('Major Loop Loss (2 Outer Legs):       %.4f W\n', 2 * PM_outer_vol_integral * Ve_c_core);
fprintf('Minor Loop Loss (2 Outer Legs):       %.4f W\n', 2 * Pm_outer_vol_integral * Ve_c_core);
fprintf('Total Loss (2 Outer Legs):          %.4f W\n', P_outer_legs_integral);
fprintf('-------------------------------------------------------------\n');
fprintf('Major Loop Loss (Center Leg):       %.4f W\n', PM_ctr_vol_integral * Ve_T_core);
fprintf('Minor Loop Loss (Center Leg):       %.4f W\n', Pm_ctr_vol_integral * Ve_T_core);
fprintf('Total Loss (Center Leg):            %.4f W\n', P_center_leg_integral);
fprintf('=============================================================\n');
fprintf('Major Loop Loss (Total):       %.4f W\n', PM_total_integral);
fprintf('Minor Loop Loss (Total):       %.4f W\n', Pm_total_integral);
fprintf('TOTAL COUPLED INDUCTOR CORE LOSS (integral):   %.4f W\n', P_total_core_integral);
fprintf('=============================================================\n');

toc; % Stop execution timer
