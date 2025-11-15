%% A E 303 - Lab 4
% Author: Parham Khodadi
% Instructor: Xiaofeng Liu

clc;clear;close all;

%% Load Data

% NACA Data
angles = [-5, 0, 5, 10, 15, 20, 20]; % Angles of attack in degrees
naca_data = struct();

naca_files = {
    'NACA_Data/01-NACA43012A_a-05.csv'
    'NACA_Data/02-NACA43012A_a+00.csv'
    'NACA_Data/03-NACA43012A_a+05.csv'
    'NACA_Data/04-NACA43012A_a+10.csv'
    'NACA_Data/05-NACA43012A_a+15.csv'
    'NACA_Data/06-NACA43012A_a+20.csv'
    'NACA_Data/07-NACA43012A_a+20.csv'
};

for i = 1:length(angles)
    % Use curly braces to extract the file name string
    data = readmatrix(naca_files{i});
    
    % Store into struct
    naca_data(i).AoA = angles(i);
    naca_data(i).x = data(:,1);
    naca_data(i).y = data(:,2);
    naca_data(i).Cp = data(:,3);
end


% Experimental Data
exp_data = struct();

exp_files = {
    'Experimental_Data/q5AoA-5.csv'
    'Experimental_Data/q5AoA0.csv'
    'Experimental_Data/q5AoA5.csv'
    'Experimental_Data/q5AoA10.csv'
    'Experimental_Data/q5AoA15.csv'
    'Experimental_Data/q5AoA20-1.csv'
    'Experimental_Data/q5AoA20-2.csv'
};

for i = 1:length(exp_files)
    data = readmatrix(exp_files{i});
    
    % Store into struct
    exp_data(i).AoA = angles(i);
    exp_data(i).Raw = data;
    
    % Optional placeholders for future computed values
    exp_data(i).Cp = []; % Will compute later
    exp_data(i).x = []; % x/c values from port locations
    exp_data(i).y = []; % y/c values from port locations
end

% Normalize with q=0 AoA=0
exp_data_0 = struct();
exp_data_0.AoA = 0;
exp_data_0.Raw = readmatrix('Experimental_Data/q0AoA0.csv');
exp_data_0.Cp = []; % Will compute later
exp_data_0.x = []; % x/c values from port locations
exp_data_0.y = []; % y/c values from port locations

%% Port Locations
% From AE_303_Lab_4_Updated_Setup.pdf, page 10

% Lower surface ports 1–16
x_lower = [0.015, 0.029, 0.055, 0.080, 0.105, 0.157, 0.207, 0.257, ...
           0.306, 0.407, 0.507, 0.608, 0.708, 0.812, 0.912, 1.000];
y_lower = [-0.0089, -0.01156, -0.0160, -0.0190, -0.0216, -0.0262, -0.0299, ...
           -0.0330, -0.0353, -0.0393, -0.0402, -0.0389, -0.0353, -0.0259, ...
           -0.0132, 0.00];

% Upper surface ports 17–32
x_upper = [0.000, 0.013, 0.025, 0.048, 0.073, 0.097, 0.150, 0.200, ...
           0.250, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900];
y_upper = [0.000, 0.0394, 0.0514, 0.0678, 0.0794, 0.0868, 0.0933, ...
           0.0927, 0.0895, 0.0857, 0.0770, 0.0656, 0.0520, 0.0370, ...
           0.0244, 0.0125];

% Full airfoil port arrays (1–32)
x_ports = [x_lower, x_upper];
y_ports = [y_lower, y_upper];

% Save for later Cp calculation
for i = 1:length(exp_data)
    exp_data(i).x = x_ports;
    exp_data(i).y = y_ports;
end

% q = 0 & AoA = 0
exp_data_0.x = x_ports;
exp_data_0.y = y_ports;

%% Conversions and Constants
inH2O_to_Psi = 0.036126; % inH2O to psi
inHg_to_Psi = 0.491154; % inHg to psi
F_to_R = 459.67;
P_amb = 30.11 * inHg_to_Psi; % psi
T_amb = 79.5 + F_to_R; % Rankine
R_air = 1716; %ft·lb/(slug·°R)

%% Compute Experimental Cp

q_inf = 5 * inH2O_to_Psi;

cp_baseline = mean(exp_data_0.Raw(9:60, 2:4001), 2);  % From q0AoA0.csv

for i = 1:length(exp_files)
    % Read the entire CSV file into a matrix
    raw_full = readmatrix(exp_files{i});
    
    % Extract the block of numeric data:
    % Rows 9 to 60 (52 ports) and Columns 2 to 4001 (4000 samples)
    data = raw_full(9:60, 2:4001);
    
    % Compute the mean pressure for each port (average over samples)
    p_avg = mean(data, 2);
    p_corrected = p_avg - mean(cp_baseline, 2); % Subtract baseline offset

    % Use Port 33 as the freestream (static) pressure reference
    p_inf = p_avg(26);
    
    % Compute Cp for ports 1 to 32
    cp = (p_corrected(1:32) - p_inf) / q_inf;
    
    % Store the computed Cp as a column vector in exp_data
    exp_data(i).Cp = cp(:);
end

%% Compute Normal (Cn) and Axial (Ca) Force Coefficients

for i = 1:length(exp_data)
    % Extract Cp for current test (nonuniform: 1–32 ports)
    cp = exp_data(i).Cp(:);
    
    % Separate upper and lower surfaces
    cp_lower = cp(1:16);  % Ports 1–16
    cp_upper = cp(17:32); % Ports 17–32

    % Define uniform grid along the chord using 161 points
    x_uniform = linspace(min(x_lower), max(x_upper), 161)';  
    % Interpolate the measured Cp data onto the uniform grid
    cp_lower_uniform = interp1(x_lower, cp_lower, x_uniform, 'linear', 'extrap');
    cp_upper_uniform = interp1(x_upper, cp_upper, x_uniform, 'linear', 'extrap');
    
    % Interpolate y-coordinates onto the uniform grid (for slope calculations)
    y_lower_uniform = interp1(x_lower, y_lower, x_uniform, 'linear', 'extrap');
    y_upper_uniform = interp1(x_upper, y_upper, x_uniform, 'linear', 'extrap');
    
    % Compute slopes on the uniform grid
    dy_dx_lower_uniform = gradient(y_lower_uniform, x_uniform);
    dy_dx_upper_uniform = gradient(y_upper_uniform, x_uniform);
    
    % ---- Cn: Normal force coefficient ----
    % (Using cp_lower - cp_upper, because the lower surface pressure is higher.)
    cn_integrand = cp_lower_uniform - cp_upper_uniform;
    cn = simpson_integration(x_uniform, cn_integrand);
    
    % ---- Ca: Axial force coefficient ----
    ca_integrand = cp_upper_uniform .* dy_dx_upper_uniform - cp_lower_uniform .* dy_dx_lower_uniform;
    ca = simpson_integration(x_uniform, ca_integrand);
    
    % Store results
    exp_data(i).Cn = cn;
    exp_data(i).Ca = ca;
end



%% Compute Lift and Drag Coefficients

for i = 1:length(exp_data)
    % Define Angle of Attack Alpha (in degrees)
    alpha = exp_data(i).AoA;
    
    % Calculate Lift Coefficient (CL) and Drag Coefficient (CD)
    CL = exp_data(i).Cn * cosd(alpha) - exp_data(i).Ca * sind(alpha);
    CD = exp_data(i).Cn * sind(alpha) + exp_data(i).Ca * cosd(alpha);
    
    % Store the results in the exp_data structure
    exp_data(i).CL = CL;
    exp_data(i).CD = CD;
end

% Display the computed Lift and Drag Coefficients
disp('Lift and Drag Coefficients for Experimental Data:');
for i = 1:length(exp_data)
    fprintf('AoA = %g°: CL = %.4f, CD = %.4f\n', exp_data(i).AoA, exp_data(i).CL, exp_data(i).CD);
end

%% Compute Pitching Moment Coefficient (Cm_LE and Cm_ac)

xac = 0.238;  % Aerodynamic center location
yac = 0.07;   % y/c location of aerodynamic center

for i = 1:length(exp_data)
    % Separate the measured Cp for upper and lower surfaces
    cp = exp_data(i).Cp(:);
    cp_lower = cp(1:16); 
    cp_upper = cp(17:32);
    
    % Define uniform grid along the chord using 161 points
    x_uniform = linspace(min(x_lower), max(x_upper), 161)';  
    % Interpolate Cp and y data onto the uniform grid
    cp_lower_uniform = interp1(x_lower, cp_lower, x_uniform, 'linear', 'extrap');
    cp_upper_uniform = interp1(x_upper, cp_upper, x_uniform, 'linear', 'extrap');
    y_lower_uniform = interp1(x_lower, y_lower, x_uniform, 'linear', 'extrap');
    y_upper_uniform = interp1(x_upper, y_upper, x_uniform, 'linear', 'extrap');
    
    % Compute slopes on the uniform grid (of the physical surfaces)
    dy_dx_lower_uniform = gradient(y_lower_uniform, x_uniform);
    dy_dx_upper_uniform = gradient(y_upper_uniform, x_uniform);
    
    % Compute moment integrals using Simpson's rule:
    % First term: (cp_upper - cp_lower) times the moment arm (x)
    I1 = simpson_integration(x_uniform, (cp_upper_uniform - cp_lower_uniform) .* x_uniform);
    
    % Second term: Upper surface slope correction
    I2 = simpson_integration(x_uniform, cp_upper_uniform .* dy_dx_upper_uniform .* y_upper_uniform);
    
    % Third term: Lower surface slope correction
    I3 = simpson_integration(x_uniform, cp_lower_uniform .* dy_dx_lower_uniform .* y_lower_uniform);
    
    % Total moment about the leading edge
    Cm_LE = I1 + I2 - I3;
    
    % Calculate Lift (CL) and Drag (CD) using the previously computed Cn and Ca
    alpha = exp_data(i).AoA;
    CL = exp_data(i).Cn * cosd(alpha) - exp_data(i).Ca * sind(alpha);
    CD = exp_data(i).Cn * sind(alpha) + exp_data(i).Ca * cosd(alpha);
    
    % Shift moment from the leading edge to the aerodynamic center
    Cm_ac = Cm_LE + CL * xac * cosd(alpha) - CD * yac * cosd(alpha) ...
                + CL * yac * sind(alpha) + CD * xac * sind(alpha);
            
    % Store the computed moments
    exp_data(i).Cm_LE = Cm_LE;
    exp_data(i).Cm_ac = Cm_ac;
end



%% Compute Aerodynamic Coefficients from NACA Data (NACA 610 43012A)

xac = 0.238;
yac = 0.07;

for i = 1:length(naca_files)
    % Read CSV data (skip header)
    data = readmatrix(naca_files{i});
    x = data(:,1); y = data(:,2); cp = data(:,3);

    % Interpolate both surfaces onto a common x grid
    x_common = linspace(0, 1, 160)';  % Match resolution of input

    % --- Correctly split surfaces from XFOIL output ---
    % XFOIL order: LE → upper surface → TE → lower surface → LE
    n = floor(length(x)/2);  % Halfway split

    x_u = x(1:n);           y_u = y(1:n);           cp_u = cp(1:n);
    x_l = flipud(x(n+1:end)); y_l = flipud(y(n+1:end)); cp_l = flipud(cp(n+1:end));

    % Interpolate Cp and y values to common grid
    cp_u_i = interp1(x_u, cp_u, x_common, 'linear', 'extrap');
    cp_l_i = interp1(x_l, cp_l, x_common, 'linear', 'extrap');
    y_u_i  = interp1(x_u, y_u, x_common, 'linear', 'extrap');
    y_l_i  = interp1(x_l, y_l, x_common, 'linear', 'extrap');

    % Compute slopes
    dy_dx_u = gradient(y_u_i) ./ gradient(x_common);
    dy_dx_l = gradient(y_l_i) ./ gradient(x_common);

    % Force coefficients
    cn = trapz(x_common, cp_l_i - cp_u_i);
    ca = trapz(x_common, cp_u_i .* dy_dx_u - cp_l_i .* dy_dx_l);

    alpha = angles(i);
    CL = cn * cosd(alpha) - ca * sind(alpha);
    CD = cn * sind(alpha) + ca * cosd(alpha);

    % Moment integrals
    I1 = trapz(x_common, (cp_u_i - cp_l_i) .* x_common);
    I2 = trapz(x_common, cp_u_i .* dy_dx_u .* y_u_i);
    I3 = trapz(x_common, cp_l_i .* dy_dx_l .* y_l_i);
    Cm_LE = I1 + I2 - I3;

    % Cm at aerodynamic center
    Cm_ac = Cm_LE + CL * xac * cosd(alpha) - CD * yac * cosd(alpha) ...
                    + CL * yac * sind(alpha) + CD * xac * sind(alpha);

    % Store to struct if needed (optional)
    naca_data(i).CL = CL;
    naca_data(i).CD = CD;
    naca_data(i).Cm_ac = Cm_ac;

    % Print results
    fprintf('AoA = %g°: CL = %.4f, CD = %.4f, Cm_ac = %.4f\n', ...
        alpha, CL, CD, Cm_ac);
end


%% Plotting Aerodynamic Coefficients vs Angle of Attack

% Experimental Data
AoA_exp = [exp_data.AoA];
CL_exp  = [exp_data.CL];
CD_exp  = [exp_data.CD];
Cm_exp  = [exp_data.Cm_ac];

% Theoretical (NACA/XFOIL) Data
AoA_theory = [naca_data.AoA];
CL_theory  = [naca_data.CL];
CD_theory  = [naca_data.CD];
Cm_theory  = [naca_data.Cm_ac];

% --- Plot CL vs AoA ---
figure;
plot(AoA_exp, CL_exp, 'o-r', 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
plot(AoA_theory, CL_theory, 's-b', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on;
xlabel('Angle of Attack (°)', 'FontSize', 12);
ylabel('Lift Coefficient (C_L)', 'FontSize', 12);
title('Coefficient of Lift vs. Angle of Attack', 'FontSize', 14);
legend('Experimental Data', 'NACA/XFOIL Data', 'Location', 'northwest');
print('-depsc2', 'CL_vs_AoA.eps');

% --- Plot CD vs AoA ---
figure;
plot(AoA_exp, CD_exp, 'o-r', 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
plot(AoA_theory, CD_theory, 's-b', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on;
xlabel('Angle of Attack (°)', 'FontSize', 12);
ylabel('Drag Coefficient (C_D)', 'FontSize', 12);
title('Coefficient of Drag vs. Angle of Attack', 'FontSize', 14);
legend('Experimental Data', 'NACA/XFOIL Data', 'Location', 'northwest');
print('-depsc2', 'CD_vs_AoA.eps');

% --- Plot Cm_ac vs AoA ---
figure;
plot(AoA_exp, Cm_exp, 'o-r', 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
plot(AoA_theory, Cm_theory, 's-b', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on;
xlabel('Angle of Attack (°)', 'FontSize', 12);
ylabel('Moment Coefficient (C_{m,ac})', 'FontSize', 12);
title('Pitching Moment Coefficient vs. Angle of Attack', 'FontSize', 14);
legend('Experimental Data', 'NACA/XFOIL Data', 'Location', 'northeast');
print('-depsc2', 'Cm_ac_vs_AoA.eps');

%% Plot Cp Distribution for All Angles of Attack

% Loop through each AoA case
figure;
hold on;
colors = lines(length(exp_data)); % distinguish curves

for i = 1:length(exp_data)
    % Extract data
    x = exp_data(i).x(:);
    cp = exp_data(i).Cp(:);

    % Plot upper surface only (ports 17–32)
    x_u = x(17:32);
    cp_u = cp(17:32);

    plot(x_u, cp_u, 'o-', 'LineWidth', 1.2, 'Color', colors(i,:), ...
         'DisplayName', sprintf('AoA = %d°', exp_data(i).AoA));
end

% Flip y-axis (by convention, lower Cp = more suction)
set(gca, 'YDir','reverse');

grid on;
xlabel('x/c', 'FontSize', 12);
ylabel('C_p (upper surface)', 'FontSize', 12);
title('Cp Distribution on Upper Surface for All AoAs', 'FontSize', 14);
legend('Location', 'best');
print('-depsc2', 'Cp_Distribution_UpperSurface.eps');


%% Compare Drag from Wake Survey vs Surface Pressure Measurement

% Constants
wake_angles = [50, 20, 15, 20, 30, 90, 90];  % degrees from AE303 Lab Setup
delta_y_prime = 0.5 * 1/12;  % 0.5 in in feet

for i = 1:length(exp_data)
    raw = exp_data(i).Raw;

    % Wake rake pressures (Ports 41–60)
    wake_block = raw(41:60, 2:end);  % 20 ports x many samples
    p_rake = mean(wake_block, 2);  % average pressure per port

    % Static & total pressures
    p_static = mean(raw(7, 2:end));
    p_total  = mean(raw(8, 2:end));

    % Dynamic pressure q1 = freestream
    q1 = p_total - p_static;

    % q2 = p_total_rake - p_static (assuming same static port for all)
    q2 = p_rake - p_static;

    % Ratio
    q2_q1 = q2 / q1;

    % Apply formula:
    % Cd = (2 / c) * ∫ (sqrt(q2/q1) - q2/q1) dy
    integrand = real(sqrt(q2_q1) - q2_q1);

    % Get vertical spacing for this AoA:
    theta_rad = deg2rad(wake_angles(i));
    dy = delta_y_prime * sin(theta_rad);  % vertical distance between ports

    % Original wake grid from 20 ports:
    y_wake = (0:19)' * dy;  % This gives 20 points (19 segments), which is not allowed by Simpson's rule

    % Create a new uniform grid with an odd number of points (e.g., 21 points)
    y_wake_new = linspace(y_wake(1), y_wake(end), 21)';
    % Interpolate the integrand onto the new grid
    integrand_new = interp1(y_wake, integrand, y_wake_new, 'linear', 'extrap');

    % Chord length in feet (assumed constant)
    c = 1.0;  % ft

    % Compute wake drag coefficient using Simpson's rule on the new grid:
    CD_wake = 2 / c * simpson_integration(y_wake_new, integrand_new);

    exp_data(i).CD_wake = CD_wake;
end



% ------------------------
% Plot comparison
% ------------------------
AoA = [exp_data.AoA];
CD_pressure = [exp_data.CD];
CD_wake = [exp_data.CD_wake];

figure;
plot(AoA, CD_pressure, 'o-r', 'LineWidth', 1.8, 'DisplayName', 'Surface Pressure');
hold on;
plot(AoA, CD_wake, 's-b', 'LineWidth', 1.8, 'DisplayName', 'Wake Survey');
grid on;
xlabel('Angle of Attack (°)', 'FontSize', 12);
ylabel('Drag Coefficient (C_D)', 'FontSize', 12);
title('Comparison of Drag Coefficients: Surface Pressure vs. Wake Survey', 'FontSize', 14);
legend('Location', 'northwest');
print('-depsc2', 'CD_Comparison_Wake_vs_Surface.eps');

%% Plot Wake Velocity Profiles for All Angles of Attack

% Wake rake info from AE_303_Lab_4_Updated_Setup.pdf
delta_y_prime = 0.5 * 1/12;  % Hypotenuse spacing between wake ports (feet)

% Angles from the table (degrees)
rake_angles = [50, 20, 15, 20, 30, 90, 90];  % One per AoA test
colors = lines(length(exp_data));

figure;
hold on;

for i = 1:length(exp_data)
    raw = exp_data(i).Raw;

    % Wake rake pressures: rows 41–60 (ports 33–52)
    p_rake = mean(raw(41:60, 2:4001), 2);  % Average over samples

    % Total pressure for test (row 8), static pressure (row 7)
    p_total = mean(raw(8, 2:4001));
    p_static = mean(raw(7, 2:4001));

    % Internal temperature (row 5 or 6 depending on your spreadsheet layout)
    T_internal = raw(5,5) + F_to_R;  % F → R

    % Density
    rho = P_amb / (R_air * T_internal);

    % Compute wake velocity (Bernoulli)
    u_wake = sqrt(2 * (p_total - p_rake) / rho);  % Vector of 20 values

    % Each port is offset by delta_y_prime along rake → vertical spacing:
    y_wake = (0:19)' * delta_y_prime * sind(rake_angles(i));  % 20 x 1 vector

    % Plot
    plot(u_wake, y_wake, 'LineWidth', 1.6, ...
        'Color', colors(i,:), ...
        'DisplayName', sprintf('AoA = %d°', exp_data(i).AoA));
end

grid on;
xlabel('Wake Velocity u (ft/s)', 'FontSize', 12);
ylabel('Vertical Position y (in)', 'FontSize', 12);
title('Wake Velocity Profiles for All AoAs', 'FontSize', 14);
legend('Location', 'best');
print('-depsc2', 'Wake_Velocity_Profiles.eps');

%% Functions
function I = simpson_integration(x, y)
    % Simpson's composite integration assumes an even number of segments.
    n = length(x);
    if mod(n-1,2) ~= 0
        error('Simpson''s rule requires an even number of segments (odd number of points).');
    end
    h = (x(end)-x(1))/(n-1);
    I = y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2));
    I = I * h / 3;
end
