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

%% PARAMETERS
lens.er = [11.9 4 2];
Ntheta = 1200;

%% COORDINATE GRID
theta_inc = linspace(eps, pi / 2, Ntheta);

%% REFRACTIVE INDECIES
lens.n = sqrt(lens.er);
lens.eccin = 1 ./ lens.n;

%% CRITICAL ANGLE AND ANGULAR DOMAIN
lens.crit_angle = asin(1 ./ lens.n);
lens.theta_max = pi / 2 - lens.crit_angle;

%% INCIDENT ANGLE @ ELIPSE
elipse.theta_inc = acos( (1 - lens.eccin' * cos(theta_inc)) ...
    ./ sqrt(1 + lens.eccin' .^ 2 - 2 * lens.eccin' * cos(theta_inc)) );
elipse.theta_inc(theta_inc > lens.theta_max') = NaN;

%% TRANSMISSION ANGLE
% @ Plane
plane.theta_tr = asin(lens.n' .* sin(theta_inc));
plane.theta_tr(imag(plane.theta_tr) ~= 0) = NaN;
% @ Elipse
elipse.theta_tr = asin(lens.n' .* sin(elipse.theta_inc));
elipse.theta_tr(imag(elipse.theta_tr) ~= 0) = NaN;

% @ Plane
plane.par_tr = NaN(length(lens.er), Ntheta);
plane.per_tr = NaN(length(lens.er), Ntheta);
plane.per_power = NaN(length(lens.er), Ntheta);
plane.par_power = NaN(length(lens.er), Ntheta);
% @ Elipse
elipse.par_tr = NaN(length(lens.er), Ntheta);
elipse.per_tr = NaN(length(lens.er), Ntheta);
elipse.per_power = NaN(length(lens.er), Ntheta);
elipse.par_power = NaN(length(lens.er), Ntheta);
for medium_idx = 1 : 1 : length(lens.er)
    %% TRANSMISSION COEFFICIENTS
    % @ Plane
    [plane.par_tr(medium_idx, :), plane.per_tr(medium_idx, :)] = ...
        transm_coeff(theta_inc, plane.theta_tr(medium_idx, :), ...
        lens.er(medium_idx), 1);
    % @ Elipse
    [elipse.par_tr(medium_idx, :), elipse.per_tr(medium_idx, :)] = ...
        transm_coeff(elipse.theta_inc(medium_idx, :), ...
        elipse.theta_tr(medium_idx, :), lens.er(medium_idx), 1);

    %% TRANSMITTED POWER
    % @ Plane
    [plane.per_power(medium_idx, :), plane.par_power(medium_idx, :)] = ...
        surf_transm_power(plane.par_tr(medium_idx, :), ...
        plane.per_tr(medium_idx, :), theta_inc, ...
        plane.theta_tr(medium_idx, :), lens.er(medium_idx), 1);
    plane.par_power(medium_idx, ...
        find(isnan(plane.par_power(medium_idx, :)), 1)) = 0;
    plane.per_power(medium_idx, ...
        find(isnan(plane.per_power(medium_idx, :)), 1)) = 0;
    % @ Elipse
    [elipse.per_power(medium_idx, :), elipse.par_power(medium_idx, :)] = ...
        surf_transm_power(elipse.par_tr(medium_idx, :), ...
        elipse.per_tr(medium_idx, :), elipse.theta_inc(medium_idx, :), ...
        elipse.theta_tr(medium_idx, :), lens.er(medium_idx), 1);
    elipse.par_power(medium_idx, ...
        find(isnan(elipse.par_power(medium_idx, :)), 1)) = 0;
    elipse.per_power(medium_idx, ...
        find(isnan(elipse.per_power(medium_idx, :)), 1)) = 0;
end

%% PLOT TRANSMISSION POWER
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4940 0.1840 0.5560];
% @ Plane
figure('Position', [250 250 1050 400]);
for medium_idx = 1 : 1 : length(lens.er)
    plot(theta_inc * 180 / pi, plane.per_power(medium_idx, :), ...
        'Color', colors(medium_idx, :), 'LineWidth', 2.0, ...
        'DisplayName',  ['TE, \epsilon_{r} = ' ...
        num2str(lens.er(medium_idx))]);
    hold on;
    plot(theta_inc * 180 / pi, plane.par_power(medium_idx, :), '--', ...
        'Color', colors(medium_idx, :), 'LineWidth', 2.0, ...
        'DisplayName', ['TM, \epsilon_{r} = ' ...
        num2str(lens.er(medium_idx))]);
    hold on;
    xline(lens.crit_angle(medium_idx) * 180 / pi, ':', 'LineWidth', 2.0, ...
        'Color', colors(medium_idx, :), 'DisplayName', ['\theta_{crit} ' ...
        '= ' num2str(round(lens.crit_angle(medium_idx) * 180 / pi, 2)) ...
        ' deg, \epsilon_{r} = ' num2str(lens.er(medium_idx))]);
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
for medium_idx = 1 : 1 : length(lens.er)
    plot(theta_inc * 180 / pi, elipse.per_power(medium_idx, :), ...
        'Color', colors(medium_idx, :), 'LineWidth', 2.0, 'DisplayName', ...
        ['TE, \epsilon_{r} = ' num2str(lens.er(medium_idx))]);
    hold on;
    plot(theta_inc * 180 / pi, elipse.par_power(medium_idx, :), '--', ...
        'Color', colors(medium_idx, :), 'LineWidth', 2.0, 'DisplayName', ...
        ['TM, \epsilon_{r} = ' num2str(lens.er(medium_idx))]);
    hold on;
    xline(lens.theta_max(medium_idx) * 180 / pi, ':', 'LineWidth', 2.0, ...
        'Color', colors(medium_idx, :), 'DisplayName', ['\theta_{max} ' ...
        '= ' num2str(round(lens.theta_max(medium_idx) * 180 / pi, 2)) ...
        ' deg, \epsilon_{r} = ' num2str(lens.er(medium_idx))]);
    hold on;
end
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('P_{t} / P_{i}');
title('TE & TM Transmitted Power Ratio @ Elipse Interface');
saveas(gcf, 'figures\tx_power_ratio_elipse.fig');

%% SAVE WORKSPACE
plane.theta_inc = theta_inc;
save('results\transmission_coeffs.mat', 'lens', 'plane', 'elipse');
