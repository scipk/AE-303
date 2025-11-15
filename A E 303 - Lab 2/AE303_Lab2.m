%% AE 303 Lab 2: Test Section Flow Uniformity Characterization
% San Diego State University
% Author: Parham Khodadi
% Date: 02/04/2025

clear; clc; close all;

%% Constants for Conversions
P_ambient_inHg = 30.16; % Ambient pressure in inHg
T_ambient_F = 71.7; % Ambient temperature in Fahrenheit
T_tunnel_F = [78.3731,76.6152,82.7341]; % Wind tunnel temperature in Fahrenheit

% Convert ambient pressure to psi (1 inHg = 0.4912 psi)
P_ambient_psi = P_ambient_inHg * 0.4912;

% Convert temperatures to Kelvin
T_ambient_K = (T_ambient_F - 32) * (5/9) + 273.15;
T_tunnel_K = (T_tunnel_F - 32) * (5/9) + 273.15;

% Conversion Factors
psi_to_Pa = 6894.76;  % Convert psi to Pascals
Pa_to_psi = 1 / psi_to_Pa;
inH2O_to_psi = 0.0360912;

%% Load Pitot Tube Geometry
geometry_file = 'PitotTube_Geometry.csv';
pitot_geometry = readmatrix(geometry_file);

% Extract X and Z coordinates for plotting
port_numbers = pitot_geometry(:, 1);
x_coords = pitot_geometry(:, 2);
z_coords = pitot_geometry(:, 3);

fprintf('Loaded Pitot Tube Geometry: %d ports.\n', length(port_numbers));

%% Compute Indicated Airspeed (IAS) for Each Test Condition

% Define test condition files
files = {'0inH2O.csv', '2inH2O.csv', '5inH2O.csv'};

% Constants
rho_air = 1.225; % Air density in kg/m^3 (standard conditions)
psi_to_Pa = 6894.76; % Conversion factor from psi to Pascals

% Initialize storage for IAS values
IAS_values = zeros(length(files), 1);
IAS_values(1) = 0;

for i = 2:length(files)
    % Read CSV file
    data = readmatrix(files{i});
    
    % Extract total and static pressures from Ports 61 & 62
    P_static_psi = data(2:601, 37);  % Static pressure at test section inlet (Port 61)
    P_total_psi = data(2:601, 38);   % Total pressure at test section inlet (Port 62)

    % Compute dynamic pressure at test section inlet (psi)
    q_psi = P_total_psi - P_static_psi;

    % Convert dynamic pressure to Pascals
    q_Pa = q_psi * psi_to_Pa;

    % Compute IAS using Bernoulli's equation: V = sqrt(2q/rho)
    airspeed_mps = sqrt(2 * q_Pa / rho_air);

    % Compute mean IAS across all 600 samples
    IAS_values(i) = mean(airspeed_mps);

    % Display results
    fprintf('Indicated Airspeed (IAS) for %s: %.3f m/s\n', files{i}, IAS_values(i));
end

% Save results for later use in analysis
save('IAS_Values.mat', 'IAS_values');


%% Objective #1 - Load Data and Compute Mean Dynamic Pressure

% Define test condition files
files = {'0inH2O.csv', '2inH2O.csv', '5inH2O.csv'};

% Initialize storage for mean dynamic pressures
mean_dynamic_pressure = zeros(length(files), 1);

% Loop through each file
for i = 1:length(files)
    % Read CSV file
    data = readmatrix(files{i});
    
    % Compute dynamic pressure: q = P_total - P_static
    dynamic_pressure_psi = data(2:601, 38) - data(2:601, 37);  % Now in psi

    % Compute mean dynamic pressure across all ports
    mean_dynamic_pressure(i) = mean(dynamic_pressure_psi(:));

    % Display results
    fprintf('Mean Dynamic Pressure for %s: %.3f psi\n', files{i}, mean_dynamic_pressure(i));
end

% Save results
save('Mean_Dynamic_Pressure.mat', 'mean_dynamic_pressure');


%% Objective #2 - Unit Reynolds Number Calculation
% This section computes the unit Reynolds number for each wind tunnel speed.

% Constants
R_specific = 287.05; % Specific gas constant for air (J/kg·K)
mu_air = 1.846e-5;   % Dynamic viscosity of air at standard conditions (Pa·s)
m_to_in = 1 / 39.3701; % Conversion factor from meters to inches

% Initialize storage for Reynolds numbers
unit_reynolds_numbers = zeros(length(files), 1);

% Compute Reynolds number for each condition
for i = 1:length(files)
    % Compute air density using the ideal gas law: rho = P / (R * T)
    rho_air = (P_ambient_psi * psi_to_Pa) / (R_specific * T_tunnel_K(i)); % in kg/m³
    
    % Compute freestream velocity: V = sqrt(2q / rho)
    V_inf = sqrt(2 * (mean_dynamic_pressure(i) * psi_to_Pa) / rho_air); % in m/s

    % Compute unit Reynolds number in 1/m
    unit_reynolds_numbers_m = (rho_air * V_inf) / mu_air; % in 1/m

    % Convert to 1/in
    unit_reynolds_numbers(i) = unit_reynolds_numbers_m * m_to_in; % in 1/in

    % Display results
    fprintf('Unit Reynolds Number for %s: Re_L = %.3e (1/in)\n', files{i}, unit_reynolds_numbers(i));
end

% Save results
save('Unit_Reynolds_Number.mat', 'unit_reynolds_numbers');

%% Objective #3 - Data Reduction & Analysis
%% Figure 1 - Averaged Pressure for Three Runs (Scatter Plot + Horizontal Averages)
figure;
hold on;

% Define test condition files
files = {'0inH2O.csv', '2inH2O.csv', '5inH2O.csv'}; 
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]}; % Blue, Orange, Yellow
legend_labels = {'$0.0$ inH$_{2}$O', '$2.0$ inH$_{2}$O', '$5.0$ inH$_{2}$O'}; % Correct LaTeX formatting

% Initialize storage for overall average pressure per dataset
all_pressures = zeros(length(files), 35); % Assuming 35 ports

for i = 1:length(files)
    % Read CSV file
    data = readmatrix(files{i});
    
    % Extract gauge pressure readings from ports 1-35 (columns 2-36)
    gauge_pressure_psi = data(2:601, 2:36); 
    
    % Compute mean pressure per port across all samples (600 samples)
    avg_pressure = mean(gauge_pressure_psi, 1); 
    
    % Store for overall averaging
    all_pressures(i, :) = avg_pressure;
    
    % Scatter plot for individual pressure measurements (no lines)
    scatter(1:35, avg_pressure, 40, colors{i}, 'filled', 'DisplayName', legend_labels{i});
    
    % Plot a horizontal black line for the average pressure of each dataset
    yline(mean(avg_pressure), 'k-', 'LineWidth', 1.5);
end

xlabel('Pressure Port', 'Interpreter', 'Latex');
ylabel('Pressure (psi)', 'Interpreter', 'Latex');
title('Flow Uniformity - Measured Pressure', 'Interpreter', 'Latex');
legend('Location', 'northeast', 'Interpreter', 'Latex');
grid on;
hold off;
print -depsc Fig1.eps


%% Figure 2 - Averaged Measured Dynamic Pressure for Three Runs
figure;
hold on;

% Define test condition files
files = {'0inH2O.csv', '2inH2O.csv', '5inH2O.csv'}; 
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]}; % Blue, Orange, Yellow
legend_labels = {'$0.0$ inH$_{2}$O', '$2.0$ inH$_{2}$O', '$5.0$ inH$_{2}$O', 'average', 'q$_{setting}$'}; % LaTeX format

% Initialize storage for dynamic pressure data
all_dynamic_pressures = zeros(length(files), 35); % Assuming 35 ports
scatter_handles = gobjects(length(files), 1); % Store scatter handles for legend

for i = 1:length(files)
    % Read CSV file
    data = readmatrix(files{i});
    
    % Extract total and static pressures **per port**
    P_total_psi = data(2:601, 2:36);  % Total pressure from ports 1-35
    P_static_psi = data(2:601, 37);   % Static pressure (single column)
    
    % Compute dynamic pressure for each port
    dynamic_pressure_psi = P_total_psi - P_static_psi; % Element-wise subtraction
    
    % Compute mean dynamic pressure per port across all samples (600 samples)
    avg_dynamic_pressure = mean(dynamic_pressure_psi, 1);
    
    % Store for overall averaging
    all_dynamic_pressures(i, :) = avg_dynamic_pressure;
    
    % Scatter plot for individual dynamic pressure measurements
    scatter(1:35, avg_dynamic_pressure, 40, colors{i}, 'o', 'filled'); % Plot points
    
    % Create a **separate scatter object for the legend** (only one per dataset)
    scatter_handles(i) = scatter(NaN, NaN, 40, colors{i}, 'o', 'filled'); % Invisible point for legend
    
    % Plot a horizontal black line for the average dynamic pressure of each dataset
    avg_line = yline(mean(avg_dynamic_pressure), 'k-', 'LineWidth', 1.5);
end

% Plot a horizontal red line for the q_setting of each dataset
q_lines = [
    yline(0, 'r-', 'LineWidth', 1.5);
    yline(2 * inH2O_to_psi, 'r-', 'LineWidth', 1.5);
    yline(5 * inH2O_to_psi, 'r-', 'LineWidth', 1.5)
];

% Create dummy scatter objects for legend entries "average" and "q_setting"
scatter_avg = scatter(NaN, NaN, 40, 'k', 's', 'filled');  % Black square for "average"
scatter_qsetting = scatter(NaN, NaN, 40, 'r', 's', 'filled');  % Red square for "q_setting"

% Set legend with scatter handles + dummy markers for average and q_setting
legend([scatter_handles; scatter_avg; scatter_qsetting], legend_labels, 'Location', 'northeast', 'Interpreter', 'Latex');

xlabel('Pressure Port', 'Interpreter', 'Latex');
ylabel('Pressure (psi)', 'Interpreter', 'Latex');
title('Flow Uniformity - Dynamic Pressure', 'Interpreter', 'Latex');
grid on;
hold off;
print -depsc Fig2.eps

%% Figure 3 - Test Section Airspeed
figure;
hold on;

load('IAS_Values.mat', 'IAS_values');

% Define test condition files
files = {'0inH2O.csv', '2inH2O.csv', '5inH2O.csv'}; 
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]}; % Blue, Orange, Yellow
legend_labels = {sprintf('%.1f m/s', IAS_values(1)), sprintf('%.1f m/s', IAS_values(2)), sprintf('%.1f m/s', IAS_values(3)), 'average', '$V_{\mathrm{setting}}$'}; % LaTeX format

rho_air = 1.225; % Air density in kg/m^3 (standard conditions)

% Initialize storage for computed airspeeds
all_airspeeds = zeros(length(files), 35); % Assuming 35 ports
scatter_handles = gobjects(length(files), 1); % Store scatter handles for legend

for i = 1:length(files)
    % Read CSV file
    data = readmatrix(files{i});
    
    % Extract total and static pressures **per port** (already in psi)
    P_total_psi = data(2:601, 2:36);  % Total pressure from ports 1-35
    P_static_psi = data(2:601, 37);   % Static pressure (single column)

    % Compute dynamic pressure for each port (psi)
    q_psi = P_total_psi - P_static_psi; % Bernoulli: q = P_total - P_static

    % Convert dynamic pressure from psi to Pascals
    q_Pa = q_psi * psi_to_Pa;

    % Fix: Ensure dynamic pressure is non-negative to avoid complex values
    q_Pa = max(q_Pa, 0); % Clamp negative values to zero

    % Compute airspeed using Bernoulli's equation: V = sqrt(2q/rho)
    airspeed_mps = sqrt(2 * q_Pa / rho_air);

    % Compute mean airspeed per port across all samples (600 samples)
    avg_airspeed = mean(airspeed_mps, 1);

    % Store for overall averaging
    all_airspeeds(i, :) = avg_airspeed;

    % Scatter plot for measured airspeed at each port
    scatter_handles(i) = scatter(1:35, avg_airspeed, 40, colors{i}, 'o', 'filled'); % Plot points
end

% Plot horizontal black line for the average airspeed of each dataset
for i = 1:length(files)
    if all(isreal(all_airspeeds(i, :))) % Ensure all values are real
        yline(mean(all_airspeeds(i, :)), 'k-', 'LineWidth', 1.5);
    end
end

% Plot red horizontal lines for indicated airspeed (IAS) settings
q_setting_lines = [
    yline(IAS_values(1), 'r-', 'LineWidth', 1.5);
    yline(IAS_values(2), 'r-', 'LineWidth', 1.5);
    yline(IAS_values(3), 'r-', 'LineWidth', 1.5)
];

% Create dummy scatter objects for legend entries "average" and "V_setting"
scatter_avg = scatter(NaN, NaN, 40, 'k', 's', 'filled');  % Black square for "average"
scatter_qsetting = scatter(NaN, NaN, 40, 'r', 's', 'filled');  % Red square for "V_setting"

% Set legend with the stored representative scatter handles
legend([scatter_handles; scatter_avg; scatter_qsetting], legend_labels, 'Location', 'northeast', 'Interpreter', 'Latex');

xlabel('Pressure Port', 'Interpreter', 'Latex');
ylabel('Airspeed (m/s)', 'Interpreter', 'Latex');
title('Flow Uniformity - Test Section Airspeed', 'Interpreter', 'Latex');
grid on;
hold off;
print -depsc Fig3.eps

%% Figure 4 - Dynamic Pressure Deviation from Average
figure;
hold on;

% Define datasets (Only 2 inH2O and 5 inH2O)
test_conditions = {'2inH2O.csv', '5inH2O.csv'};
colors = {[0, 0, 0], [1, 0, 0]}; % Black for 2 inH2O, Red for 5 inH2O
legend_labels = {'$2.0$ inH$_{2}$O', '$5.0$ inH$_{2}$O'}; % LaTeX format

% Load mean dynamic pressure values from Objective 1
load('Mean_Dynamic_Pressure.mat', 'mean_dynamic_pressure');

% Define port numbers (excluding Port 7)
ports = 1:35; 
ports(7) = []; % Remove Port 7 (static pressure probe)

% Initialize storage for deviations
dq_values = zeros(length(test_conditions), length(ports));

for i = 1:length(test_conditions)
    % Read CSV file
    data = readmatrix(test_conditions{i});
    
    % Extract total and static pressures for ports 1-35
    P_total_psi = data(2:601, 2:36);  % Total pressure (ports 1-35)
    P_static_psi = data(2:601, 37);   % Static pressure (single column)

    % Compute dynamic pressure: q = P_total - P_static
    dynamic_pressure_psi = P_total_psi - P_static_psi;

    % Compute mean dynamic pressure per port
    avg_dynamic_pressure = mean(dynamic_pressure_psi, 1);

    % Compute deviation from **mean dynamic pressure per port**
    deviation = 100 * (avg_dynamic_pressure - mean_dynamic_pressure(i+1)) ./ mean_dynamic_pressure(i+1);

    deviation(7) = [];
    % Store for plotting
    dq_values(i, :) = deviation;
end

% Ensure correct data dimensions
if size(dq_values, 2) ~= length(ports)
    error('Mismatch in vector sizes! Check ports and computed deviations.');
end

% Plot deviations
scatter(ports, dq_values(1, :), 50, 'MarkerEdgeColor', colors{1}, 'Marker', 'o', 'MarkerFaceColor', 'none'); % Black for 2 inH2O
scatter(ports, dq_values(2, :), 50, 'MarkerEdgeColor', colors{2}, 'Marker', 'd', 'MarkerFaceColor', 'none'); % Red for 5 inH2O

% Set legend
legend(legend_labels, 'Location', 'southeast', 'Interpreter', 'Latex');

% Set labels and title
xlabel('Pressure Port Number', 'Interpreter', 'Latex');
ylabel('$\frac{q - \bar{q}}{\bar{q}} \times 100\%$', 'Interpreter', 'Latex');
title('Flow Uniformity - Dynamic Pressure Deviation From Average', 'Interpreter', 'Latex');

% Formatting
ytickformat('percentage');
xlim([1 35]);
ylim([-5 2.5]); % Match expected range
xticks(1:2:35);
yticks(-5:0.5:2.5);
grid on;
hold off;
print -depsc Fig4.eps

%% Figure 5 - Dynamic Pressure Deviation for 2 inH2O
figure;
hold on;

% Load mean dynamic pressure values from Objective 1
load('Mean_Dynamic_Pressure.mat', 'mean_dynamic_pressure');

% Define dataset and properties
test_condition = '2inH2O.csv';
color = [0, 0, 0]; % Black
marker = 'o'; % Circle
legend_label = '$2.0$ inH$_{2}$O'; % LaTeX format

% Define port numbers (excluding Port 7)
ports = 1:35; 
ports(7) = []; % Remove Port 7 (static pressure probe)

% Read CSV file
data = readmatrix(test_condition);

% Extract total and static pressures
P_total_psi = data(2:601, 2:36);  % Total pressure (ports 1-35)
P_static_psi = data(2:601, 37);   % Static pressure (single column)

% Compute dynamic pressure: q = P_total - P_static
dynamic_pressure_psi = P_total_psi - P_static_psi;

% Compute mean dynamic pressure per port
avg_dynamic_pressure = mean(dynamic_pressure_psi, 1);

% Compute deviation from **mean dynamic pressure per port**
deviation = 100 * (avg_dynamic_pressure - mean_dynamic_pressure(2)) ./ mean_dynamic_pressure(2);

% Remove Port 7 data
deviation(7) = [];

% Scatter plot
scatter(ports, deviation, 50, 'MarkerEdgeColor', color, 'Marker', marker, 'MarkerFaceColor', 'none');

% Set legend
legend(legend_label, 'Location', 'southeast', 'Interpreter', 'Latex');

% Set labels and title
xlabel('Pressure Port Number', 'Interpreter', 'Latex');
ylabel('$\frac{q - \bar{q}}{\bar{q}} \times 100\%$', 'Interpreter', 'Latex');
title('Flow Uniformity - Dynamic Pressure Deviation from Average 2.0 inH$_2$O', 'Interpreter', 'Latex');

% Formatting
ytickformat('percentage');
xlim([1 35]);
ylim([-5 2.5]); % Match expected range
xticks(1:2:35);
yticks(-5:0.5:2.5);
grid on;
hold off;
print -depsc Fig5.eps

%% Figure 6 - Dynamic Pressure Deviation for 5 inH2O
figure;
hold on;

% Define dataset and properties
test_condition = '5inH2O.csv';
color = [1, 0, 0]; % Red
marker = 'd'; % Diamond
legend_label = '$5.0$ inH$_{2}$O'; % LaTeX format

% Read CSV file
data = readmatrix(test_condition);

% Extract total and static pressures
P_total_psi = data(2:601, 2:36);  % Total pressure (ports 1-35)
P_static_psi = data(2:601, 37);   % Static pressure (single column)

% Compute dynamic pressure: q = P_total - P_static
dynamic_pressure_psi = P_total_psi - P_static_psi;

% Compute mean dynamic pressure per port
avg_dynamic_pressure = mean(dynamic_pressure_psi, 1);

% Compute deviation from **mean dynamic pressure per port**
deviation = 100 * (avg_dynamic_pressure - mean_dynamic_pressure(3)) ./ mean_dynamic_pressure(3);

% Remove Port 7 data
deviation(7) = [];

% Scatter plot
scatter(ports, deviation, 50, 'MarkerEdgeColor', color, 'Marker', marker, 'MarkerFaceColor', 'none');

% Set legend
legend(legend_label, 'Location', 'southeast', 'Interpreter', 'Latex');

% Set labels and title
xlabel('Pressure Port Number', 'Interpreter', 'Latex');
ylabel('$\frac{q - \bar{q}}{\bar{q}} \times 100\%$', 'Interpreter', 'Latex');
title('Flow Uniformity - Dynamic Pressure Deviation from Average 5.0 inH$_2$O', 'Interpreter', 'Latex');

% Formatting
ytickformat('percentage');
xlim([1 35]);
ylim([-5 2.5]); % Match expected range
xticks(1:2:35);
yticks(-5:0.5:2.5);
grid on;
hold off;
print -depsc Fig6.eps


%% Figures 7-9 - Compute and Save Mean Airspeed

% Define test condition files
files = {'0inH2O.csv', '2inH2O.csv', '5inH2O.csv'};

% Constants
rho_air = 1.225; % Air density in kg/m^3 (standard conditions)
psi_to_Pa = 6894.76; % Convert psi to Pascals

% Initialize storage for mean airspeed
mean_airspeed = zeros(length(files), 1);

% Loop through each test condition
for i = 1:length(files)
    % Read CSV file
    data = readmatrix(files{i});
    
    % Compute dynamic pressure: q = P_total - P_static
    dynamic_pressure_psi = data(2:601, 38) - data(2:601, 37); % q in psi

    % Convert dynamic pressure to Pascals
    q_Pa = dynamic_pressure_psi * psi_to_Pa;

    % Ensure q is non-negative to prevent complex numbers
    q_Pa = max(q_Pa, 0); 

    % Compute airspeed using Bernoulli's equation: V = sqrt(2q / rho)
    airspeed_mps = sqrt(2 * q_Pa / rho_air);

    % Compute mean airspeed across all ports
    mean_airspeed(i) = mean(airspeed_mps(:));

    % Display results
    fprintf('Mean Airspeed for %s: %.3f m/s\n', files{i}, mean_airspeed(i));
end

% Save results for later use
save('Mean_Airspeed.mat', 'mean_airspeed');

%% Figure 7 - Airspeed Deviation for 2 inH2O and 5 inH2O
figure;
hold on;

% Define test conditions (2 inH2O and 5 inH2O)
test_conditions = {'2inH2O.csv', '5inH2O.csv'};
colors = {[0, 0, 0], [1, 0, 0]}; % Black for 2 inH2O, Red for 5 inH2O
markers = {'o', 'd'}; % Circle for 2 inH2O, Diamond for 5 inH2O
legend_labels = {'$2.0$ inH$_{2}$O', '$5.0$ inH$_{2}$O'}; % LaTeX format

% Load mean airspeed values
load('Mean_Airspeed.mat', 'mean_airspeed');

% Define port numbers (excluding Port 7)
ports = 1:35; 
ports(7) = []; % Remove Port 7 (static pressure probe)

% Initialize storage for airspeed deviations
V_deviation = zeros(length(test_conditions), length(ports));

for i = 1:length(test_conditions)
    % Read CSV file
    data = readmatrix(test_conditions{i});
    
    % Extract total and static pressures for ports 1-35
    P_total_psi = data(2:601, 2:36);  % Total pressure (ports 1-35)
    P_static_psi = data(2:601, 37);   % Static pressure (single column)

    % Compute dynamic pressure: q = P_total - P_static
    dynamic_pressure_psi = P_total_psi - P_static_psi;

    % Convert dynamic pressure to Pascals
    q_Pa = dynamic_pressure_psi * psi_to_Pa;

    % Ensure q is non-negative
    q_Pa = max(q_Pa, 0); 

    % Compute airspeed using Bernoulli’s equation
    airspeed = sqrt(2 * q_Pa / rho_air);

    % Compute mean airspeed per port
    avg_airspeed = mean(airspeed, 1);

    % Compute deviation from **mean airspeed per port**
    deviation = 100 * (avg_airspeed - mean_airspeed(i+1)) ./ mean_airspeed(i+1);

    % Remove Port 7 data
    deviation(7) = [];

    % Store for plotting
    V_deviation(i, :) = deviation;
end

% Scatter plot for 2 inH2O and 5 inH2O
scatter(ports, V_deviation(1, :), 50, 'MarkerEdgeColor', colors{1}, 'Marker', 'o', 'MarkerFaceColor', 'none'); % Black for 2 inH2O
scatter(ports, V_deviation(2, :), 50, 'MarkerEdgeColor', colors{2}, 'Marker', 'd', 'MarkerFaceColor', 'none'); % Red for 5 inH2O

% Set legend
legend(legend_labels, 'Location', 'southeast', 'Interpreter', 'Latex');

% Set labels and title
xlabel('Pressure Port Number', 'Interpreter', 'Latex');
ylabel('$\frac{V - \bar{V}}{\bar{V}} \times 100$ (\%)', 'Interpreter', 'Latex');
title('Flow Uniformity - Airspeed Deviation from Average', 'Interpreter', 'Latex');

% Formatting
ytickformat('percentage');
xlim([1 35]);
ylim([-2.5 1.5]); % Match expected range
xticks(1:2:35);
yticks(-2.5:0.5:1.5);
grid on;
hold off;
print -depsc Fig7.eps

%% Figure 8 - Airspeed Deviation for 2 inH2O
figure;
hold on;

% Scatter plot for 2 inH2O
scatter(ports, V_deviation(1, :), 50, 'MarkerEdgeColor', colors{1}, 'Marker', 'o', 'MarkerFaceColor', 'none');

% Set legend
legend(legend_labels{1}, 'Location', 'southeast', 'Interpreter', 'Latex');

% Set labels and title
xlabel('Pressure Port Number', 'Interpreter', 'Latex');
ylabel('$\frac{V - \bar{V}}{\bar{V}} \times 100$ (\%)', 'Interpreter', 'Latex');
title('Flow Uniformity - Airspeed Deviation from Average 2.0 inH$_2$O', 'Interpreter', 'Latex');

% Formatting
ytickformat('percentage');
xlim([1 35]);
ylim([-2.5 1.5]); % Match expected range
xticks(1:2:35);
yticks(-2.5:0.5:1.5);
grid on;
hold off;
print -depsc Fig8.eps

%% Figure 9 - Airspeed Deviation for 5 inH2O
figure;
hold on;

% Scatter plot for 5 inH2O
scatter(ports, V_deviation(2, :), 50, 'MarkerEdgeColor', colors{2}, 'Marker', 'd', 'MarkerFaceColor', 'none');

% Set legend
legend(legend_labels{2}, 'Location', 'southeast', 'Interpreter', 'Latex');

% Set labels and title
xlabel('Pressure Port Number', 'Interpreter', 'Latex');
ylabel('$\frac{V - \bar{V}}{\bar{V}} \times 100$ (\%)', 'Interpreter', 'Latex');
title('Flow Uniformity - Airspeed Deviation from Average 5.0 inH$_2$O', 'Interpreter', 'Latex');

% Formatting
ytickformat('percentage');
xlim([1 35]);
ylim([-2.5 1.5]); % Match expected range
xticks(1:2:35);
yticks(-2.5:0.5:1.5);
grid on;
hold off;
print -depsc Fig9.eps

%% Figures 10 & 11 - Dynamic Pressure Deviation Contour Maps

% Load Pitot Tube Geometry
geometry_file = 'PitotTube_Geometry.csv';
pitot_geometry = readmatrix(geometry_file);
port_numbers = pitot_geometry(:, 1);
x_coords = pitot_geometry(:, 2);
z_coords = pitot_geometry(:, 3);

% Remove Port 7 from dataset
valid_ports = port_numbers ~= 7;
x_coords = x_coords(valid_ports);
z_coords = z_coords(valid_ports);
port_numbers = port_numbers(valid_ports);

% Define test condition files (2 inH2O and 5 inH2O)
files = {'2inH2O.csv', '5inH2O.csv'};
titles = {...
    'Dynamic pressure deviation from average $q=2.0$ inH$_2$O', ...
    'Dynamic pressure deviation from average $q=5.0$ inH$_2$O' ...
};
colorbar_labels = '$\frac{q - \bar{q}}{\bar{q}} \times 100\%$';

% Load mean dynamic pressure from Objective #1
load('Mean_Dynamic_Pressure.mat', 'mean_dynamic_pressure');

% Create a grid for interpolation
xq = linspace(min(x_coords), max(x_coords), 50);
zq = linspace(min(z_coords), max(z_coords), 50);
[Xq, Zq] = meshgrid(xq, zq);

for fig = 1:2
    % Read pressure data
    data = readmatrix(files{fig});

    % Extract total and static pressures for ports 1-35
    P_total_psi = data(2:601, 2:36);  % Total pressure from ports 1-35
    P_static_psi = data(2:601, 37);   % Static pressure (single column)

    % Compute dynamic pressure per port
    dynamic_pressure_psi = P_total_psi - P_static_psi;

    % Compute mean dynamic pressure per port
    avg_dynamic_pressure = mean(dynamic_pressure_psi, 1);

    % Compute dynamic pressure deviation from mean
    deviation = 100 * (avg_dynamic_pressure - mean_dynamic_pressure(fig+1)) ./ mean_dynamic_pressure(fig+1);

    % Remove Port 7 from deviation data
    deviation = deviation(valid_ports);

    % Interpolate data onto grid for smooth contour plot
    deviation_interp = griddata(x_coords, z_coords, deviation, Xq, Zq, 'cubic');

    % Create new figure window for each test condition
    figure;
    hold on;
    
    % Contour plot with black stroke between levels
    contourf(Xq, Zq, deviation_interp, 20, 'LineColor', 'k'); % Add black strokes
    scatter(x_coords, z_coords, 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); % Black-edged Pitot locations
    
    % Formatting
    colorbar_handle = colorbar;
    title(titles{fig}, 'Interpreter', 'Latex', 'FontSize', 14);
    xlabel('x position [in]', 'Interpreter', 'Latex', 'FontSize', 12);
    ylabel('z position [in]', 'Interpreter', 'Latex', 'FontSize', 12);
    
    % Correct LaTeX formatting for colorbar label
    colorbar_handle.Label.String = colorbar_labels;
    colorbar_handle.Label.Interpreter = 'Latex';
    colorbar_handle.Label.FontSize = 12;
    
    hold off;
    print(sprintf('Fig%i.eps', 9+fig),'-depsc')
end


%% Figure 12 - Airspeed Contour Map for Highest IAS Setting

% Load Pitot Tube Geometry
geometry_file = 'PitotTube_Geometry.csv';
pitot_geometry = readmatrix(geometry_file);
port_numbers = pitot_geometry(:, 1);
x_coords = pitot_geometry(:, 2);
z_coords = pitot_geometry(:, 3);

% Remove Port 7 from dataset
valid_ports = port_numbers ~= 7;
x_coords = x_coords(valid_ports);
z_coords = z_coords(valid_ports);
port_numbers = port_numbers(valid_ports);

% Load Indicated Airspeed (IAS) values
load('IAS_Values.mat', 'IAS_values');
IAS_target = IAS_values(3); % Third IAS value (highest setting)

% Load the 5inH2O dataset (corresponding to highest IAS)
file = '5inH2O.csv';
data = readmatrix(file);

% Extract total and static pressures for ports 1-35
P_total_psi = data(2:601, 2:36);  % Total pressure from ports 1-35
P_static_psi = data(2:601, 37);   % Static pressure (single column)

% Compute dynamic pressure per port
dynamic_pressure_psi = P_total_psi - P_static_psi;

% Constants
rho_air = 1.225; % Air density in kg/m^3 (standard conditions)
psi_to_Pa = 6894.76; % Conversion factor from psi to Pascals

% Convert dynamic pressure to Pascals
q_Pa = dynamic_pressure_psi * psi_to_Pa;

% Compute airspeed using Bernoulli's equation: V = sqrt(2q/rho)
airspeed_mps = sqrt(2 * q_Pa / rho_air);

% Compute mean airspeed per port
avg_airspeed = mean(airspeed_mps, 1);

% Remove Port 7 from airspeed data
avg_airspeed = avg_airspeed(valid_ports);

% Create a grid for interpolation
xq = linspace(min(x_coords), max(x_coords), 50);
zq = linspace(min(z_coords), max(z_coords), 50);
[Xq, Zq] = meshgrid(xq, zq);

% Interpolate airspeed data onto grid
airspeed_interp = griddata(x_coords, z_coords, avg_airspeed, Xq, Zq, 'cubic');

% Create figure
figure;
hold on;

% Contour plot with black stroke between levels
contourf(Xq, Zq, airspeed_interp, 20, 'LineColor', 'k'); % Add black strokes
scatter(x_coords, z_coords, 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); % Black-edged Pitot locations

% Formatting
colorbar_handle = colorbar;
title(sprintf('Airspeed distribution - IAS=%.1f m/s', IAS_target), 'Interpreter', 'Latex', 'FontSize', 14);
xlabel('x position [in]', 'Interpreter', 'Latex', 'FontSize', 12);
ylabel('z position [in]', 'Interpreter', 'Latex', 'FontSize', 12);

% Correct LaTeX formatting for colorbar label
colorbar_handle.Label.String = 'Measured Airspeed (m/s)';
colorbar_handle.Label.Interpreter = 'Latex';
colorbar_handle.Label.FontSize = 12;

hold off;
print -depsc Fig12.eps

%% Figure 13 - Airspeed Contour Map for 2 inH2O Setting

% Load Pitot Tube Geometry
geometry_file = 'PitotTube_Geometry.csv';
pitot_geometry = readmatrix(geometry_file);
port_numbers = pitot_geometry(:, 1);
x_coords = pitot_geometry(:, 2);
z_coords = pitot_geometry(:, 3);

% Remove Port 7 from dataset
valid_ports = port_numbers ~= 7;
x_coords = x_coords(valid_ports);
z_coords = z_coords(valid_ports);
port_numbers = port_numbers(valid_ports);

% Load Indicated Airspeed (IAS) values
load('IAS_Values.mat', 'IAS_values');
IAS_target = IAS_values(2); % Second IAS value (corresponding to 2 inH2O)

% Load the 2inH2O dataset
file = '2inH2O.csv';
data = readmatrix(file);

% Extract total and static pressures for ports 1-35
P_total_psi = data(2:601, 2:36);  % Total pressure from ports 1-35
P_static_psi = data(2:601, 37);   % Static pressure (single column)

% Compute dynamic pressure per port
dynamic_pressure_psi = P_total_psi - P_static_psi;

% Constants
rho_air = 1.225; % Air density in kg/m^3 (standard conditions)
psi_to_Pa = 6894.76; % Conversion factor from psi to Pascals

% Convert dynamic pressure to Pascals
q_Pa = dynamic_pressure_psi * psi_to_Pa;

% Compute airspeed using Bernoulli's equation: V = sqrt(2q/rho)
airspeed_mps = sqrt(2 * q_Pa / rho_air);

% Compute mean airspeed per port
avg_airspeed = mean(airspeed_mps, 1);

% Remove Port 7 from airspeed data
avg_airspeed = avg_airspeed(valid_ports);

% Create a grid for interpolation
xq = linspace(min(x_coords), max(x_coords), 50);
zq = linspace(min(z_coords), max(z_coords), 50);
[Xq, Zq] = meshgrid(xq, zq);

% Interpolate airspeed data onto grid
airspeed_interp = griddata(x_coords, z_coords, avg_airspeed, Xq, Zq, 'cubic');

% Create figure
figure;
hold on;

% Contour plot with black stroke between levels
contourf(Xq, Zq, airspeed_interp, 20, 'LineColor', 'k'); % Black strokes for levels
scatter(x_coords, z_coords, 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); % Black-edged Pitot locations

% Formatting
colorbar_handle = colorbar;
title(sprintf('Airspeed distribution - IAS=%.1f m/s', IAS_target), 'Interpreter', 'Latex', 'FontSize', 14);
xlabel('x position [in]', 'Interpreter', 'Latex', 'FontSize', 12);
ylabel('z position [in]', 'Interpreter', 'Latex', 'FontSize', 12);

% Correct LaTeX formatting for colorbar label
colorbar_handle.Label.String = 'Measured Airspeed (m/s)';
colorbar_handle.Label.Interpreter = 'Latex';
colorbar_handle.Label.FontSize = 12;

hold off;
print -depsc Fig13.eps


%% Objective #4 - Convergence Test for Time-Averaged Dynamic Pressure

% Define test condition files
files = {'2inH2O.csv', '5inH2O.csv'};
colors = {'b', 'r'}; % Blue for 2inH2O, Red for 5inH2O
labels = {'2 inH_2O', '5 inH_2O'};
threshold = 0.005 / 100; % Convert 0.005% to decimal

figure;
hold on;

for i = 1:length(files)
    % Read CSV file
    data = readmatrix(files{i});
    
    % Compute dynamic pressure (q = P_total - P_static)
    q_psi = data(2:601, 38) - data(2:601, 37); % Dynamic pressure from all 600 samples
    
    % Compute cumulative mean at each sample index
    cumulative_mean = cumsum(q_psi) ./ (1:length(q_psi))';

    % Compute convergence criteria: |q_avg(n) - q_avg(n-1)| / q_avg(n-1)
    convergence_error = abs(diff(cumulative_mean) ./ cumulative_mean(2:end)) * 100; % Convert to %

    % Find the **last** sample index where error exceeds threshold
    last_convergence_index = find(convergence_error > threshold * 100, 1, 'last') + 1;

    % If no valid index is found, assume full sample range
    if isempty(last_convergence_index)
        last_convergence_index = length(q_psi);
    end

    % Plot convergence curve
    plot(2:length(q_psi), convergence_error, 'Color', colors{i}, 'LineWidth', 1.2, 'DisplayName', labels{i});
    
    % Plot horizontal threshold line
    yline(threshold * 100, '--k', 'LineWidth', 1.5, 'Label', '0.005%', 'Interpreter', 'Latex');
    
    % Plot vertical marker for the last convergence sample
    xline(last_convergence_index, '--', 'Color', colors{i}, 'LineWidth', 1.5, 'Label', sprintf('Last: %d', last_convergence_index));
    
    % Display result
    fprintf('Dynamic pressure converges for %s at %d samples.\n', files{i}, last_convergence_index);
end

% Formatting
xlabel('Number of Samples', 'Interpreter', 'Latex');
ylabel('Convergence Error (\%)', 'Interpreter', 'Latex');
title('Dynamic Pressure Convergence', 'Interpreter', 'Latex');
legend('show', 'Location', 'northeast', 'Interpreter', 'Latex');
ytickformat('percentage'); % Format y-axis in %
grid on;
hold off;
print -depsc Fig14.eps
