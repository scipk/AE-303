%% AE 303 Lab 6
% Author: Parham Khodadi
% Instructor: Xiaofeng Liu

clear; clc; close all;

%% Load Data
% Set working directory to data folder
cd('EduPIV_Lab_Data');

% Run the import script
run('ImportData_EduPIV.m');

% Return to parent folder for analysis output
cd('..');

%% Data Processing: Time-Averaged Fields, Reynolds Stresses, and Vorticity

% Nozzle diameter in mm (given as 5 cm)
D = 50;

% Time-averaged velocities (raw is more complete, avoids NaNs)
u_mean = mean(u_raw, 3, 'omitnan');
v_mean = mean(v_raw, 3, 'omitnan');

% Fluctuations (Reynolds decomposition from filtered velocities)
u_fluct = bsxfun(@minus, u, mean(u, 3, 'omitnan'));
v_fluct = bsxfun(@minus, v, mean(v, 3, 'omitnan'));

% Reynolds stresses from filtered vectors
uu = mean(u_fluct.^2, 3, 'omitnan');
vv = mean(v_fluct.^2, 3, 'omitnan');
uv = mean(u_fluct .* v_fluct, 3, 'omitnan');

% --- Uniform Grid Generation ---
x_vec = linspace(min(X(:)), max(X(:)), size(X,2));
y_vec = linspace(min(Y(:)), max(Y(:)), size(Y,1));
[X_uniform, Y_uniform] = meshgrid(x_vec, y_vec);

% Interpolated Reynolds stresses
uu_uniform = griddata(X, Y, uu, X_uniform, Y_uniform, 'linear');
vv_uniform = griddata(X, Y, vv, X_uniform, Y_uniform, 'linear');
uv_uniform = griddata(X, Y, uv, X_uniform, Y_uniform, 'linear');

% Interpolate u_mean and v_mean onto the uniform grid
u_uniform = griddata(X, Y, u_mean, X_uniform, Y_uniform, 'linear');
v_uniform = griddata(X, Y, v_mean, X_uniform, Y_uniform, 'linear');

% Fill any remaining NaNs at the edges
u_uniform = fillmissing(u_uniform, 'nearest');
v_uniform = fillmissing(v_uniform, 'nearest');

% Compute gradients and vorticity
dx = mean(diff(x_vec));
dy = mean(diff(y_vec));
[du_dy, ~] = gradient(u_uniform, dy, dx);
[~, dv_dx] = gradient(v_uniform, dy, dx);
vorticity = dv_dx - du_dy;

% Non-dimensional coordinates
X_nd = X_uniform / D;
Y_nd = Y_uniform / D;

%% Objective 2: Time-Averaged u-Velocity Contour Map

figure;
contourf(X_nd, Y_nd, u_uniform, 20, 'LineStyle', 'none');
axis equal;
colorbar;
xlabel('x/D');
ylabel('y/D');
title('Time-Averaged u Velocity Contour Map (\itū\rm)');
saveas(gcf, 'figs/u_mean_contour.eps', 'epsc'); % Save .eps file

%% Objective 3: Reynolds Normal Stress Contour Maps (\itu'^2\rm and \itv'^2\rm)

% u'2 contour
figure;
contourf(X_nd, Y_nd, uu_uniform, 20, 'LineStyle', 'none');
axis equal;
colorbar;
xlabel('x/D'); ylabel('y/D');
title('Reynolds Normal Stress Contour: \itu''^2\rm');
saveas(gcf, 'figs/uu_reynolds_normal_stress.eps', 'epsc'); % Save .eps file

% v'2 contour
figure;
contourf(X_nd, Y_nd, vv_uniform, 20, 'LineStyle', 'none');
axis equal;
colorbar;
xlabel('x/D'); ylabel('y/D');
title('Reynolds Normal Stress Contour: \itv''^2\rm');
saveas(gcf, 'figs/vv_reynolds_normal_stress.eps', 'epsc'); 

%% Objective 4: Reynolds Shear Stress Contour Map (\itu'v'\rm)

figure;
contourf(X_nd, Y_nd, uv_uniform, 20, 'LineStyle', 'none');
axis equal;
colorbar;
xlabel('x/D'); ylabel('y/D');
title('Reynolds Shear Stress Contour: \itu''v''\rm');

% Save .eps file
saveas(gcf, 'figs/uv_reynolds_shear_stress.eps', 'epsc');

%% Objective 5: Vorticity Magnitude Contour Map

% Compute magnitude of vorticity
vort_mag = abs(vorticity);

figure;
contourf(X_nd, Y_nd, vort_mag, 20, 'LineStyle', 'none');
axis equal;
colorbar;
xlabel('x/D'); ylabel('y/D');
title('Vorticity Magnitude Contour');

% Save .eps file
saveas(gcf, 'figs/vorticity_magnitude.eps', 'epsc');

%% Objective 6: Profiles at Selected Streamwise Locations

% Define five x/D positions to extract vertical profiles
xD_locs = [0.05, 0.10, 0.15, 0.20, 0.25];

% Initialize color options for visual distinction
colors = lines(length(xD_locs));

% Function to find the closest column index for each x/D location
get_x_index = @(x_val) find(abs(X_nd(1,:) - x_val) == min(abs(X_nd(1,:) - x_val)), 1);

% Plot: Mean u-velocity profiles
figure; hold on;
for i = 1:length(xD_locs)
    col = get_x_index(xD_locs(i));
    plot(u_uniform(:,col), Y_nd(:,col), 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('x/D = %.2f', xD_locs(i)));
end
xlabel('u (m/s)');
ylabel('y/D');
title('Mean u Velocity Profiles at Selected x/D Locations');
legend;
grid on;
saveas(gcf, 'figs/u_profiles_vs_y.eps', 'epsc');

% Plot: Reynolds normal stress u''²
figure; hold on;
for i = 1:length(xD_locs)
    col = get_x_index(xD_locs(i));
    plot(uu_uniform(:,col), Y_nd(:,col), 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('x/D = %.2f', xD_locs(i)));
end
xlabel('u''^2 (m^2/s^2)');
ylabel('y/D');
title('Reynolds Normal Stress (u''^2) Profiles');
legend;
grid on;
saveas(gcf, 'figs/uu_profiles_vs_y.eps', 'epsc');

% Plot: Reynolds normal stress v''²
figure; hold on;
for i = 1:length(xD_locs)
    col = get_x_index(xD_locs(i));
    plot(vv_uniform(:,col), Y_nd(:,col), 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('x/D = %.2f', xD_locs(i)));
end
xlabel('v''^2 (m^2/s^2)');
ylabel('y/D');
title('Reynolds Normal Stress (v''^2) Profiles');
legend;
grid on;
saveas(gcf, 'figs/vv_profiles_vs_y.eps', 'epsc');

% Plot: Reynolds shear stress u''v''
figure; hold on;
for i = 1:length(xD_locs)
    col = get_x_index(xD_locs(i));
    plot(uv_uniform(:,col), Y_nd(:,col), 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('x/D = %.2f', xD_locs(i)));
end
xlabel('u''v'' (m^2/s^2)');
ylabel('y/D');
title('Reynolds Shear Stress Profiles');
legend;
grid on;
saveas(gcf, 'figs/uv_profiles_vs_y.eps', 'epsc');

% Plot: Vorticity profiles
figure; hold on;
for i = 1:length(xD_locs)
    col = get_x_index(xD_locs(i));
    plot(vorticity(:,col), Y_nd(:,col), 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', sprintf('x/D = %.2f', xD_locs(i)));
end
xlabel('\omega_z (1/s)');
ylabel('y/D');
title('Vorticity Profiles at Selected x/D Locations');
legend;
grid on;
saveas(gcf, 'figs/vorticity_profiles_vs_y.eps', 'epsc');
