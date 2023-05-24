close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

if ~exist([pwd() '\results'], 'dir')
    mkdir('results');
end

addpath('../quasi-optics-library');

lens(3) = struct('er', [], 'n', [], 'e', [], 'theta_crit', [], ...
    'theta_max', []);
ellipse(3) = struct('theta_inc', [], 'theta_tr', [], 'TM_coef', [], ...
    'TE_coef', [], 'TM_ratio', [], 'TE_ratio', []);
plane(3) = struct('theta_inc', [], 'theta_tr', [], 'TM_coef', [], ...
    'TE_coef', [], 'TM_ratio', [], 'TE_ratio', []);

%% PARAMETERS
% Lens parameters
lens(1).er = 11.9;
lens(2).er = 4;
lens(3).er = 2;
% Grid parameters
Ntheta = 1200;

%% DEPENDENT PARAMETERS
% Lens parameters
for lens_idx = 1 : 1 : length(lens)
    lens(lens_idx).n = sqrt(lens(lens_idx).er);
    lens(lens_idx).e = 1 ./ lens(lens_idx).n;
    lens(lens_idx).theta_crit = asin(1 ./ lens(lens_idx).n);
    lens(lens_idx).theta_max = pi / 2 - lens(lens_idx).theta_crit;
end

%% COORDINATE GRID
theta_inc = linspace(eps, pi / 2, Ntheta);

%% INCIDENT ANGLE @ ELIPSE
for lens_idx = 1 : 1 : length(lens)
    ellipse(lens_idx).theta_inc = acos( (1 - lens(lens_idx).e * cos(theta_inc)) ...
        ./ sqrt(1 + lens(lens_idx).e .^ 2 - 2 * lens(lens_idx).e * cos(theta_inc)) );
    ellipse(lens_idx).theta_inc(theta_inc > lens(lens_idx).theta_max') = NaN;
end

%% TRANSMISSION ANGLE
for lens_idx = 1 : 1 : length(lens)
    % @ Plane
    plane(lens_idx).theta_tr = asin(lens(lens_idx).n .* sin(theta_inc));
    plane(lens_idx).theta_tr(imag(plane(lens_idx).theta_tr) ~= 0) = NaN;
    % @ Elipse
    ellipse(lens_idx).theta_tr = asin(lens(lens_idx).n ...
        .* sin(ellipse(lens_idx).theta_inc));
    ellipse(lens_idx).theta_tr(imag(ellipse(lens_idx).theta_tr) ~= 0) = NaN;
end

for lens_idx = 1 : 1 : length(lens)
    %% TRANSMISSION COEFFICIENTS
    % @ Plane
    [plane(lens_idx).TM_coef, plane(lens_idx).TE_coef] = ...
        transm_coeff(theta_inc, plane(lens_idx).theta_tr, ...
        lens(lens_idx).er, 1);
    % @ Elipse
    [ellipse(lens_idx).TM_coef, ellipse(lens_idx).TE_coef] = ...
        transm_coeff(ellipse(lens_idx).theta_inc, ...
        ellipse(lens_idx).theta_tr, lens(lens_idx).er, 1);

    %% TRANSMITTED POWER
    % @ Plane
    [plane(lens_idx).TE_ratio, plane(lens_idx).TM_ratio] = ...
        surf_transm_power(plane(lens_idx).TM_coef, ...
        plane(lens_idx).TE_coef, theta_inc, ...
        plane(lens_idx).theta_tr, lens(lens_idx).er, 1);
    plane(lens_idx).TM_ratio(find(isnan(plane(lens_idx).TM_ratio), 1)) = 0;
    plane(lens_idx).TE_ratio(find(isnan(plane(lens_idx).TE_ratio), 1)) = 0;
    % @ Elipse
    [ellipse(lens_idx).TE_ratio, ellipse(lens_idx).TM_ratio] = ...
        surf_transm_power(ellipse(lens_idx).TM_coef, ...
        ellipse(lens_idx).TE_coef, ellipse(lens_idx).theta_inc, ...
        ellipse(lens_idx).theta_tr, lens(lens_idx).er, 1);
    ellipse(lens_idx).TM_ratio(find(isnan(ellipse(lens_idx).TM_ratio), 1)) = 0;
    ellipse(lens_idx).TE_ratio(find(isnan(ellipse(lens_idx).TE_ratio), 1)) = 0;
end

%% PLOT TRANSMISSION POWER
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4940 0.1840 0.5560];
% @ Plane
figure('Position', [250 250 1050 400]);
for lens_idx = 1 : 1 : length(lens)
    plot(theta_inc * 180 / pi, plane(lens_idx).TE_ratio, ...
        'Color', colors(lens_idx, :), 'LineWidth', 2.0, ...
        'DisplayName',  ['TE, \epsilon_{r} = ' ...
        num2str(lens(lens_idx).er)]);
    hold on;
    plot(theta_inc * 180 / pi, plane(lens_idx).TM_ratio, '--', ...
        'Color', colors(lens_idx, :), 'LineWidth', 2.0, ...
        'DisplayName', ['TM, \epsilon_{r} = ' ...
        num2str(lens(lens_idx).er)]);
    hold on;
    xline(lens(lens_idx).theta_crit * 180 / pi, ':', 'LineWidth', 2.0, ...
        'Color', colors(lens_idx, :), 'DisplayName', ['\theta_{crit} ' ...
        '= ' num2str(round(lens(lens_idx).theta_crit * 180 / pi, 2)) ...
        ' deg, \epsilon_{r} = ' num2str(lens(lens_idx).er)]);
    hold on;
end
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('P_{t} / P_{i}');
title('TE & TM Transmitted Power Ratio @ Plane Interface');
saveas(gcf, 'figures\tx_power_ratio_plane.fig');
% @ Elipse
figure('Position', [250 250 1050 400]);
for lens_idx = 1 : 1 : length(lens)
    plot(theta_inc * 180 / pi, ellipse(lens_idx).TE_ratio, ...
        'Color', colors(lens_idx, :), 'LineWidth', 2.0, 'DisplayName', ...
        ['TE, \epsilon_{r} = ' num2str(lens(lens_idx).er)]);
    hold on;
    plot(theta_inc * 180 / pi, ellipse(lens_idx).TM_ratio, '--', ...
        'Color', colors(lens_idx, :), 'LineWidth', 2.0, 'DisplayName', ...
        ['TM, \epsilon_{r} = ' num2str(lens(lens_idx).er)]);
    hold on;
    xline(lens(lens_idx).theta_max * 180 / pi, ':', 'LineWidth', 2.0, ...
        'Color', colors(lens_idx, :), 'DisplayName', ['\theta_{max} ' ...
        '= ' num2str(round(lens(lens_idx).theta_max * 180 / pi, 2)) ...
        ' deg, \epsilon_{r} = ' num2str(lens(lens_idx).er)]);
    hold on;
end
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('P_{t} / P_{i}');
title('TE & TM Transmitted Power Ratio @ Ellipse Interface');
saveas(gcf, 'figures\tx_power_ratio_ellipse.fig');

%% SAVE WORKSPACE
for plane_idx = 1 : 1 : length(plane)
    plane(plane_idx).theta_inc = theta_inc;
end
save('results\transmission_coeffs.mat', 'lens', 'plane', 'ellipse');
