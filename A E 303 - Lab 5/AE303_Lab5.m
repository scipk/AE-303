%% AE 303 Lab 5: Full Model Aircraft Aerodynamic Analysis
% Author: Parham Khodadi
% Instructor: Xiaofeng Liu

clear; clc; close all;

%% Load Data
data = readtable('data.csv');

%% Constants and Conversions
inH2OtoPsi = 27.7076;   % Convert psi to inH2O
psitoPa = 6894.76;      % Convert psi to Pascals
FtoR = 459.67;          % Convert from Fahrenheit to Rankine

S = 93.81;              % Reference area (in^2)
c_bar = 3.466;              % Chord length (in)
b = 27.066;             % Span (in)
q = 7 / inH2OtoPsi;     % Dynamic pressure in psi
R = 1716;               % ft*lbf / slug*R

%% Separate Runs
run1 = data(2:15,:);    % Tail off (if used)
run2 = data(17:31,:);   % Main test (Tail on)
run3 = data(33:47,:);   % Gravity tare
run4 = data(49:53,:);   % Aero tare
run5 = data(55:59,:);   % Baseline

%% Define columns
% Columns in the data
col_alpha = 5;  % Alpha
col_beta  = 6;  % Beta
col_T = 7;      % T (Temperature)
col_Fx = 8;     % Fx (Drag Force)
col_Fy = 9;     % Fy (Side Force)
col_Fz = 10;    % Fz (Lift Force)
col_Mx = 11;    % Mx (Rolling Moment)
col_My = 12;    % My (Pitching Moment)
col_Mz = 13;    % Mz (Yawing Moment)

%% Data Correction (Run 2)

run2_corrected = run2;

% Create the correction data for when there is no Model (run4 - run3)
% but only for the Fx,Fy,Fz,Mx,My,Mz columns
modeloff_correction = run4;
for i = 3:5
    modeloff_correction(i,col_Fx:col_Mz) = modeloff_correction(i,col_Fx:col_Mz) - run5(i,col_Fx:col_Mz);
end

% Correct run2 by run2_corrected = run2 - run1 - (run4 - run3)
% but only for the Fx,Fy,Fz,Mx,My,Mz columns
for i = 4:15
    if run2{i,col_beta} == 0
        run2_corrected(i,col_Fx:col_Mz) = run2_corrected(i,col_Fx:col_Mz) - run1(i-1,col_Fx:col_Mz) - modeloff_correction(3,col_Fx:col_Mz);
    elseif run2{i,col_beta} == 5
        run2_corrected(i,col_Fx:col_Mz) = run2_corrected(i,col_Fx:col_Mz) - run1(i-1,col_Fx:col_Mz) - modeloff_correction(4,col_Fx:col_Mz);
    elseif run2{i,col_beta} == 10
        run2_corrected(i,col_Fx:col_Mz) = run2_corrected(i,col_Fx:col_Mz) - run1(i-1,col_Fx:col_Mz) - modeloff_correction(5,col_Fx:col_Mz);
    end
end

%% Data Correction (Run 3)
run3_corrected = run3;

% Correct run3 by run3_corrected = run3 - (run4 - run3)
% but only for the Fx,Fy,Fz,Mx,My,Mz columns
% Using the same (run4 - run3) data from run2_corrected
for i = 4:15
    if run3{i,col_beta} == 0
        run3_corrected(i,col_Fx:col_Mz) = run3_corrected(i,col_Fx:col_Mz) - modeloff_correction(3,col_Fx:col_Mz);
    elseif run3{i,col_beta} == 5
        run3_corrected(i,col_Fx:col_Mz) = run3_corrected(i,col_Fx:col_Mz) - modeloff_correction(4,col_Fx:col_Mz);
    elseif run3{i,col_beta} == 10
        run3_corrected(i,col_Fx:col_Mz) = run3_corrected(i,col_Fx:col_Mz) - modeloff_correction(5,col_Fx:col_Mz);
    end
end

%% Calculate Aerodynamic Force and Moment Coefficients

% Angles of Attack
alpha_run2 = run2_corrected{4:15,col_alpha};

% Betas
beta_run2 = run2_corrected{4:15, col_beta};

% Lift Coefficient = L/(q*S) = Fz/(q*S)
C_L = run2_corrected{4:15,col_Fz}/(q*S);

% Drag Coefficient = D/(q*S) = Fx/(q*S)
C_D = run2_corrected{4:15, col_Fx}/(q*S);

% Pitching Moment Coefficient = M/(q*S*c_bar) = My/(q*S*c_bar)
C_M = run2_corrected{4:15, col_My}/(q*S*c_bar);

% Yaw Moment Coefficient = N/(q*S*b) = Mz/(q*S*b)
C_N = run2_corrected{4:15, col_Mz}/(q*S*b);

%% Plot Aerodynamic Coefficients

mask = beta_run2 == 0;  % for symmetric alpha sweep plots

% CL vs alpha
figure;
plot(alpha_run2(mask), C_L(mask), 'o-', 'LineWidth', 1.5);
xlabel('\alpha (deg)');
ylabel('C_L');
title('C_L vs \alpha (Tail-On)');
grid on;
print(gcf, '-depsc2', 'CL_vs_alpha.eps');

% CM vs alpha
figure;
plot(alpha_run2(mask), C_M(mask), 's-', 'LineWidth', 1.5);
xlabel('\alpha (deg)');
ylabel('C_M');
title('C_M vs \alpha (Tail-On)');
grid on;
print(gcf, '-depsc2', 'CM_vs_alpha.eps');

% CN vs beta (only use rows with alpha = 0 to isolate beta sweep)
mask_beta = run2_corrected{:, col_alpha} == 0;
beta_tailon = run2_corrected{mask_beta, col_beta};
CN_tailon = run2_corrected{mask_beta, col_Mz} / (q * S * b);

figure;
plot(beta_tailon, CN_tailon, 'd-', 'LineWidth', 1.5);
xlabel('\beta (deg)');
ylabel('C_N');
title('C_N vs \beta (Tail-On)');
grid on;
print(gcf, '-depsc2', 'CN_vs_beta.eps');

% CL vs CD
figure;
plot(C_D(mask), C_L(mask), '^-', 'LineWidth', 1.5);
xlabel('C_D');
ylabel('C_L');
title('C_L vs C_D (Tail-On)');
grid on;
print(gcf, '-depsc2', 'CL_vs_CD.eps');

%% Calculate Oswald Efficiency Factor (e)

% Aspect Ratio (from lab docs)
AR = (27.066)^2 / 93.81;

% Filter region for parabolic drag polar fit
CL_fit = C_L(mask);
CD_fit = C_D(mask);

% Only use CL in range [0.2, 0.8]
fit_mask = (CL_fit >= 0.3) & (CL_fit <= 0.6);
CL_range = CL_fit(fit_mask);
CD_range = CD_fit(fit_mask);

% Prepare design matrix: CD = CD0 + K*CL^2
X = [ones(length(CL_range),1), CL_range.^2];
coeffs = X \ CD_range;  % Linear least squares

CD0 = coeffs(1);
K = coeffs(2);

% Oswald efficiency factor
e = 1 / (pi * K * AR);

% Display results
fprintf('Parabolic Drag Polar Fit Results:\n');
fprintf('CD0 = %.4f\n', CD0);
fprintf('K = %.4f\n', K);
fprintf('Oswald Efficiency Factor e = %.4f\n', e);

% Plot the fit
CL_plot = linspace(min(CL_range), max(CL_range), 100);
CD_plot = CD0 + K * CL_plot.^2;

figure;
plot(CL_range, CD_range, 'bo', 'DisplayName', 'Data');
hold on;
plot(CL_plot, CD_plot, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Fit');
xlabel('C_L');
ylabel('C_D');
title('Drag Polar Fit: C_D = C_{D,0} + K C_L^2');
legend('Location','northwest');
grid on;
print(gcf, '-depsc2', 'DragPolarFit_CD_vs_CL.eps');

%% Find Maximum Lift-to-Drag Ratio (C_L/C_D)

% Use only data where beta = 0 (mask already defined earlier)
CL_use = C_L(mask);
CD_use = C_D(mask);

% Compute L/D for all points
L_over_D = CL_use ./ CD_use;

% Find the maximum L/D and its index
[LoverD_max, idx_max] = max(L_over_D);

% Extract values
CL_max_LD = CL_use(idx_max);
CD_max_LD = CD_use(idx_max);
alpha_max_LD = alpha_run2(mask);
alpha_at_max_LD = alpha_max_LD(idx_max);

% Display results
fprintf('\nMaximum C_L/C_D:\n');
fprintf('Max L/D = %.3f at alpha = %.2f deg\n', LoverD_max, alpha_at_max_LD);
fprintf('C_L = %.4f, C_D = %.4f\n', CL_max_LD, CD_max_LD);

% Plot L/D vs alpha
figure;
plot(alpha_max_LD, L_over_D, 'ko-', 'LineWidth', 1.5);
xlabel('\alpha (deg)');
ylabel('C_L / C_D');
title('Lift-to-Drag Ratio vs \alpha (Tail-On)');
grid on;
print(gcf, '-depsc2', 'CL_over_CD_vs_alpha.eps');


%% Find Zero-Lift Angle of Attack and Lift Curve Slope

% Filter to linear region of CL vs alpha
alpha_lin = alpha_run2(mask);
CL_lin = C_L(mask);

% Choose linear range manually (usually -6 to +6 deg)
linear_mask = (alpha_lin >= -6) & (alpha_lin <= 6);
alpha_fit = alpha_lin(linear_mask);
CL_fit = CL_lin(linear_mask);

% Linear fit: CL = a * alpha + b
coeffs_CL = polyfit(alpha_fit, CL_fit, 1);
dCL_dalpha = coeffs_CL(1);    % Slope
alpha_L0 = -coeffs_CL(2)/coeffs_CL(1);   % Intercept = -b/a

% Display results
fprintf('\nLift Curve Fit Results:\n');
fprintf('dCL/dalpha = %.4f per deg\n', dCL_dalpha);
fprintf('Zero-lift angle of attack alpha_L=0 = %.2f deg\n', alpha_L0);

% Plot
alpha_plot = linspace(-6, 6, 100);  % restrict to linear region only
CL_plot = polyval(coeffs_CL, alpha_plot);

figure;
plot(alpha_lin, CL_lin, 'bo', 'DisplayName', 'Data');
hold on;
plot(alpha_plot, CL_plot, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit');
xlabel('\alpha (deg)');
ylabel('C_L');
title('Lift Curve: Linear Fit to Find dC_L/d\alpha and \alpha_{L=0}');
legend('Location','northwest');
grid on;
print(gcf, '-depsc2', 'CL_vs_alpha_linear_fit.eps');

%% Find Maximum Lift Coefficient and Stall Angle

% Use clean alpha/C_L sweep from earlier (beta = 0)
alpha_trim = alpha_run2(mask);
CL_trim = C_L(mask);

% Find maximum CL
[CL_max, idx_max] = max(CL_trim);
alpha_stall = alpha_trim(idx_max);

% Display results
fprintf('\nMaximum Lift and Stall:\n');
fprintf('C_L,max = %.4f at alpha = %.2f deg (stall angle)\n', CL_max, alpha_stall);

% Plot
figure;
plot(alpha_trim, CL_trim, 'bo-', 'LineWidth', 1.5);
hold on;
plot(alpha_stall, CL_max, 'r*', 'MarkerSize', 10, 'DisplayName', 'C_{L,max}');
xlabel('\alpha (deg)');
ylabel('C_L');
title('Maximum Lift Coefficient and Stall Angle');
grid on;
legend('show');
print(gcf, '-depsc2', 'CLmax_vs_alpha.eps');


%% Fit dCM/dalpha for Tail-On and Tail-Off

% --- Tail-On ---
alpha_tailon = alpha_run2(mask);
CM_tailon = C_M(mask);

% Use only linear range, e.g. [-6 to 6 deg]
fit_mask_tailon = (alpha_tailon >= -6) & (alpha_tailon <= 6);
alpha_fit_tailon = alpha_tailon(fit_mask_tailon);
CM_fit_tailon = CM_tailon(fit_mask_tailon);

coeffs_CM_tailon = polyfit(alpha_fit_tailon, CM_fit_tailon, 1);
dCM_dalpha_tailon = coeffs_CM_tailon(1);

% Define alpha and C_M for tail-off configuration
alpha_run3 = run3_corrected{4:15, col_alpha};      % Angle of attack (Run 3)
C_M_tailoff = run3_corrected{4:15, col_My} / (q * S * c_bar);  % Pitching moment coefficient (Run 3)

% --- Tail-Off ---
alpha_tailoff = alpha_run3(mask);  % run3 corrected alpha
CM_tailoff = C_M_tailoff(mask);    % run3 corrected moment

fit_mask_tailoff = (alpha_tailoff >= -6) & (alpha_tailoff <= 6);
alpha_fit_tailoff = alpha_tailoff(fit_mask_tailoff);
CM_fit_tailoff = CM_tailoff(fit_mask_tailoff);

coeffs_CM_tailoff = polyfit(alpha_fit_tailoff, CM_fit_tailoff, 1);
dCM_dalpha_tailoff = coeffs_CM_tailoff(1);

% Display results
fprintf('\nPitching Moment Slope Results:\n');
fprintf('Tail-On:    dCM/dalpha = %.4f\n', dCM_dalpha_tailon);
fprintf('Tail-Off:   dCM/dalpha = %.4f\n', dCM_dalpha_tailoff);

% Plot for Tail-On
figure;
plot(alpha_tailon, CM_tailon, 'bs', 'DisplayName', 'Tail-On Data');
hold on;
plot(alpha_fit_tailon, polyval(coeffs_CM_tailon, alpha_fit_tailon), 'b-', 'LineWidth', 1.5);

% Plot for Tail-Off
plot(alpha_tailoff, CM_tailoff, 'rs', 'DisplayName', 'Tail-Off Data');
plot(alpha_fit_tailoff, polyval(coeffs_CM_tailoff, alpha_fit_tailoff), 'r-', 'LineWidth', 1.5);

xlabel('\alpha (deg)');
ylabel('C_M');
title('Pitching Moment vs \alpha: Tail-On vs Tail-Off');
legend('Location','northeast');
grid on;
print(gcf, '-depsc2', 'CM_vs_alpha_tailon_tailoff.eps');

%% Find dCN/dBeta for Tail-On and Tail-Off

% --- Tail-On ---
mask_beta_tailon = run2_corrected{:, col_alpha} == 0;  % constant alpha = 0
beta_tailon = run2_corrected{mask_beta_tailon, col_beta};
CN_tailon = run2_corrected{mask_beta_tailon, col_Mz} / (q * S * b);

% Linear fit
coeffs_CN_tailon = polyfit(beta_tailon, CN_tailon, 1);
dCN_dbeta_tailon = coeffs_CN_tailon(1);

% --- Tail-Off ---
mask_beta_tailoff = run3_corrected{:, col_alpha} == 0;
beta_tailoff = run3_corrected{mask_beta_tailoff, col_beta};
CN_tailoff = run3_corrected{mask_beta_tailoff, col_Mz} / (q * S * b);

% Linear fit
coeffs_CN_tailoff = polyfit(beta_tailoff, CN_tailoff, 1);
dCN_dbeta_tailoff = coeffs_CN_tailoff(1);

% Display results
fprintf('\nYawing Moment Slope Results:\n');
fprintf('Tail-On:    dCN/dbeta = %.4f\n', dCN_dbeta_tailon);
fprintf('Tail-Off:   dCN/dbeta = %.4f\n', dCN_dbeta_tailoff);

% Plot
figure;
plot(beta_tailon, CN_tailon, 'bd', 'DisplayName', 'Tail-On Data');
hold on;
plot(beta_tailoff, CN_tailoff, 'rd', 'DisplayName', 'Tail-Off Data');

plot(beta_tailon, polyval(coeffs_CN_tailon, beta_tailon), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Tail-On Fit');
plot(beta_tailoff, polyval(coeffs_CN_tailoff, beta_tailoff), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Tail-Off Fit');

xlabel('\beta (deg)');
ylabel('C_N');
title('Yawing Moment vs \beta: Tail-On vs Tail-Off');
legend('Location','northwest');
grid on;
print(gcf, '-depsc2', 'CN_vs_beta_tailon_tailoff.eps');
