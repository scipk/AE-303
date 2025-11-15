%% A E 303 - Lab 3 - Wind tunnel free stream turbulence measurement using a turbulence sphere
% Name: Parham Khodadi
% Instructor: Xiaofeng Liu
% San Diego State University

clc; clear; close all;

% Define file names
file_4in = '4in_Sphere.csv';
file_4987in = '4.987in_Sphere.csv';
file_6in = '6in_Sphere.csv';

% Load data
data_4in = readmatrix(file_4in);
data_4987in = readmatrix(file_4987in);
data_6in = readmatrix(file_6in);

% Display size of data to confirm import
disp('Data Loaded.');

%% Conversion Factors
inH2O_to_psi = 1/27.708;
inHg_to_psi = 1/2.036;
F_to_R = 459.67;
psi_to_lbf_ft2 = 144;
K_to_R = 1.8;


%% Load some variables

% q_setting (psi)
q_setting_6in = data_6in(6:4:62, 2) * inH2O_to_psi;
q_setting_4in = data_4in(6:4:62, 2) * inH2O_to_psi;
q_setting_4987in = data_4987in(6:4:62, 2) * inH2O_to_psi;

% T_setting (°R)
T_setting_6in = data_6in(6:4:62, 2) + F_to_R;
T_setting_4in = data_4in(6:4:62, 2) + F_to_R;
T_setting_4987in = data_4987in(6:4:62, 2) + F_to_R;

% P_ambient (psi)
P_ambient_6in = data_6in(2,3) * inHg_to_psi;
P_ambient_4in = data_4in(2,3) * inHg_to_psi;
P_ambient_4987in = data_4987in(2,3) * inHg_to_psi;

% T_ambient (°R)
T_ambient_6in = data_6in(1,3) + F_to_R;
T_ambient_4in = data_4in(1,3) + F_to_R;
T_ambient_4987in = data_4987in(1,3) + F_to_R;

%% Calculate Dynamic Pressure (q)
% Extract pressure data
P_total_4in = data_4in(7:4:63, 5:804); % Port 2: Total pressure
P_static_4in = data_4in(6:4:62, 5:804); % Port 1: Static pressure

P_total_4987in = data_4987in(7:4:63, 5:804); % Port 2: Total pressure
P_static_4987in = data_4987in(6:4:62, 5:804); % Port 1: Static pressure

P_total_6in = data_6in(7:4:63, 5:804); % Port 2: Total pressure
P_static_6in = data_6in(6:4:62, 5:804); % Port 1: Static pressure

% Compute dynamic pressure (q)
q_4in = P_total_4in - P_static_4in;
q_4987in = P_total_4987in - P_static_4987in;
q_6in = P_total_6in - P_static_6in;

disp('Calculated Dynamic Pressure.');

%% Calculate Pressure Difference (Δp)

% Extract pressure data
P_stagnation_4in = data_4in(9:4:65, 5:804); % Port 4: Stagnation pressure
P_aft_4in = data_4in(8:4:64, 5:804); % Port 3: Aft (rear) pressure

P_stagnation_4987in = data_4987in(9:4:65, 5:804); % Port 4: Stagnation pressure
P_aft_4987in = data_4987in(8:4:64, 5:804); % Port 3: Aft (rear) pressure

P_stagnation_6in = data_6in(9:4:65, 5:804); % Port 4: Stagnation pressure
P_aft_6in = data_6in(8:4:64, 5:804); % Port 3: Aft (rear) pressure

% Compute pressure difference (ΔP)
delta_P_4in = P_stagnation_4in - P_aft_4in;
delta_P_4987in = P_stagnation_4987in - P_aft_4987in;
delta_P_6in = P_stagnation_6in - P_aft_6in;

disp('Calculated Pressure Difference (ΔP).');

%% Compute Normalized Pressure Difference (ΔP/q)

% Calculate ΔP/q for each sphere
delta_P_over_q_4in = delta_P_4in ./ q_4in;
delta_P_over_q_4987in = delta_P_4987in ./ q_4987in;
delta_P_over_q_6in = delta_P_6in ./ q_6in;

disp('Calculated Normalized Pressure Difference (ΔP/q).');

%% Compute Mean ΔP/q Across 800 Samples for Each Test Point

% Compute mean for each Reynolds test point
delta_P_over_q_4in_mean = mean(delta_P_over_q_4in, 2); 
delta_P_over_q_4987in_mean = mean(delta_P_over_q_4987in, 2);
delta_P_over_q_6in_mean = mean(delta_P_over_q_6in, 2);

disp('Computed Mean ΔP/q for Each Test Point.');

%% Compute Air Density (ρ_test) using Ideal Gas Law

% Convert P_ambient from psi to lb/ft²
P_ambient_4in_lbf_ft2 = P_ambient_4in * psi_to_lbf_ft2;
P_ambient_4987in_lbf_ft2 = P_ambient_4987in * psi_to_lbf_ft2;
P_ambient_6in_lbf_ft2 = P_ambient_6in * psi_to_lbf_ft2;

% Define specific gas constant for air in ft·lb/(slug·°R)
R_air = 1716; 

% Compute air density (slug/ft³)
rho_test_4in = P_ambient_4in_lbf_ft2 ./ (R_air * T_ambient_4in);
rho_test_4987in = P_ambient_4987in_lbf_ft2 ./ (R_air * T_ambient_4987in);
rho_test_6in = P_ambient_6in_lbf_ft2 ./ (R_air * T_ambient_6in);

disp('Corrected Air Density (ρ) Calculation.');


%% Compute Air Viscosity (μ_test) using Sutherland's Formula

% Sutherland's constants
mu_0 = 1.716e-5 * 0.020885; % Reference viscosity in slug/(ft·s)
T_0 = 273 * K_to_R; % Reference temperature in Rankine
S = 111 * K_to_R; % Sutherland's constant for air in Rankine

% Compute viscosity for each sphere
mu_test_4in = mu_0 * ((T_ambient_4in / T_0).^1.5) .* ((T_0 + S) ./ (T_ambient_4in + S));
mu_test_4987in = mu_0 * ((T_ambient_4987in / T_0).^1.5) .* ((T_0 + S) ./ (T_ambient_4987in + S));
mu_test_6in = mu_0 * ((T_ambient_6in / T_0).^1.5) .* ((T_0 + S) ./ (T_ambient_6in + S));

disp('Corrected Air Viscosity (μ_test) Calculation.');

%% Compute Reynolds Number (Re_test)

% Sphere diameters in feet
D_4in = 4 / 12;
D_4987in = 4.987 / 12;
D_6in = 6 / 12;

% Compute velocity (U) using dynamic pressure q
U_4in = sqrt(2 * q_setting_4in * psi_to_lbf_ft2/ rho_test_4in);
U_4987in = sqrt(2 * q_setting_4987in * psi_to_lbf_ft2/ rho_test_4987in);
U_6in = sqrt(2 * q_setting_6in * psi_to_lbf_ft2/ rho_test_6in);

% Compute Reynolds number (Re_tunnel)
Re_test_4in = (rho_test_4in .* U_4in .* D_4in) / mu_test_4in;
Re_test_4987in = (rho_test_4987in .* U_4987in .* D_4987in) / mu_test_4987in;
Re_test_6in = (rho_test_6in .* U_6in .* D_6in) / mu_test_6in;

disp('Calculated Reynolds Number (Re_tunnel).');

%% Identify and Remove Outlier for 4-inch Sphere (11th q_setting)
outlier_index = 11; % 11th dynamic pressure setting

% Remove outlier from ΔP/q and Re_tunnel for the 4-inch sphere
Re_test_4in(outlier_index) = [];
delta_P_over_q_4in_mean(outlier_index) = [];

disp('Outlier removed from 4-inch sphere data.');


%% Re-Plot Re_tunnel vs. Mean ΔP/q Without Outlier
figure;
hold on;
grid on;
box on;

% Custom colors
color_4in = [0 0.4470 0.7410];  % Blue
color_4987in = [0.8500 0.3250 0.0980]; % Red-Orange
color_6in = [0.9290 0.6940 0.1250]; % Yellow-Gold

% Plot after removing outlier
plot(Re_test_4in, delta_P_over_q_4in_mean, 'o-', 'Color', color_4in, ...
     'MarkerFaceColor', color_4in, 'MarkerEdgeColor', color_4in, ...
     'DisplayName', 'D = 4.000 in.');
plot(Re_test_4987in, delta_P_over_q_4987in_mean, 's-', 'Color', color_4987in, ...
     'MarkerFaceColor', color_4987in, 'MarkerEdgeColor', color_4987in, ...
     'DisplayName', 'D = 4.987 in.');
plot(Re_test_6in, delta_P_over_q_6in_mean, 'd-', 'Color', color_6in, ...
     'MarkerFaceColor', color_6in, 'MarkerEdgeColor', color_6in, ...
     'DisplayName', 'D = 6.000 in.');

% Reference Line for ΔP/q = 1.220
yline(1.220, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Re_c');

% Labels and Title
xlabel('Reynolds number');
ylabel('\Deltap/q');
title('Reynolds vs \Deltap/q (Outlier Removed)');
legend('Location', 'Best');

hold off;

% Save Figure as EPS
saveas(gcf, 'Reynolds_vs_DeltaPq.eps', 'epsc');


disp('Plotted Re_tunnel vs. Mean ΔP/q without the outlier.');

%% Identify Critical Reynolds Number (Re_c)

% Define threshold for ΔP/q
threshold = 1.220;

% Find index closest to ΔP/q = 1.220 for each sphere
[~, idx_4in] = min(abs(delta_P_over_q_4in_mean - threshold));
[~, idx_4987in] = min(abs(delta_P_over_q_4987in_mean - threshold));
[~, idx_6in] = min(abs(delta_P_over_q_6in_mean - threshold));

% Extract critical Reynolds numbers
Re_c_4in = Re_test_4in(idx_4in);
Re_c_4987in = Re_test_4987in(idx_4987in);
Re_c_6in = Re_test_6in(idx_6in);

disp(['Critical Reynolds Number (Re_c) for 4 in Sphere: ', num2str(Re_c_4in)]);
disp(['Critical Reynolds Number (Re_c) for 4.987 in Sphere: ', num2str(Re_c_4987in)]);
disp(['Critical Reynolds Number (Re_c) for 6 in Sphere: ', num2str(Re_c_6in)]);

%% Compute Turbulence Factor (TF)
TF_4in = 385000 / Re_c_4in;
TF_4987in = 385000 / Re_c_4987in;
TF_6in = 385000 / Re_c_6in;

disp(['Turbulence Factor (TF) for 4 in Sphere: ', num2str(TF_4in)]);
disp(['Turbulence Factor (TF) for 4.987 in Sphere: ', num2str(TF_4987in)]);
disp(['Turbulence Factor (TF) for 6 in Sphere: ', num2str(TF_6in)]);

%% Compute Effective Reynolds Number (Re_effective)
Re_effective_4in = TF_4in * Re_test_4in;
Re_effective_4987in = TF_4987in * Re_test_4987in;
Re_effective_6in = TF_6in * Re_test_6in;

disp('Computed Effective Reynolds Number (Re_effective).');

%% Compute Average Turbulence Intensity
I_4in = 2.15; % Eyeballed from Figure 2
I_4987in = 2.05;
I_6in = 2.10;

I_avg = mean([I_4in, I_4987in, I_6in]);

% Display final results
disp(['Turbulence Intensity for 4 in Sphere: ', num2str(I_4in), '%']);
disp(['Turbulence Intensity for 4.987 in Sphere: ', num2str(I_4987in), '%']);
disp(['Turbulence Intensity for 6 in Sphere: ', num2str(I_6in), '%']);
disp(['Average Turbulence Intensity: ', num2str(I_avg), '%']);