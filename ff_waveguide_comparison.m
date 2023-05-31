close all;
clear;
clc;

ff_waveguide_silicon;
ff_waveguide_quartz;
ff_waveguide_plastic;

close all;

load('results\wg_field_silicon.mat');
load('results\wg_field_quartz.mat');
load('results\wg_field_plastic.mat');

Es_total = norm_magnitude(wg_field_silicon.E_total, 'dB');
Eq_total = norm_magnitude(wg_field_quatz.E_total, 'dB');
Ep_total = norm_magnitude(wg_field_plastic.E_total, 'dB');

theta_s_plot = NaN(1, 2 * length(wg_field_silicon.theta));
theta_s_plot(1 : length(wg_field_silicon.theta)) = - fliplr(wg_field_silicon.theta) * 180 / pi;
theta_s_plot(length(wg_field_silicon.theta) + 1 : end) = wg_field_silicon.theta * 180 / pi;
theta_q_plot = NaN(1, 2 * length(wg_field_quatz.theta));
theta_q_plot(1 : length(wg_field_quatz.theta)) = - fliplr(wg_field_quatz.theta) * 180 / pi;
theta_q_plot(length(wg_field_quatz.theta) + 1 : end) = wg_field_quatz.theta * 180 / pi;
theta_p_plot = NaN(1, 2 * length(wg_field_plastic.theta));
theta_p_plot(1 : length(wg_field_plastic.theta)) = - fliplr(wg_field_plastic.theta) * 180 / pi;
theta_p_plot(length(wg_field_plastic.theta) + 1 : end) = wg_field_plastic.theta * 180 / pi;

figure('Position', [250 250 725 400]);
plane_field = NaN(1, length(theta_s_plot));
plane_field(1 : length(wg_field_silicon.theta)) = fliplr(Es_total(1601, :));
plane_field(length(wg_field_silicon.theta) + 1 : end) = Es_total(1, :);
plot(theta_s_plot, plane_field, 'LineWidth', 2.0, ...
    'DisplayName', 'E_{WG}^{FF}, silicon');
hold on;
plane_field = NaN(1, length(theta_q_plot));
plane_field(1 : length(wg_field_quatz.theta)) = fliplr(Eq_total(1601, :));
plane_field(length(wg_field_quatz.theta) + 1 : end) = Eq_total(1, :);
plot(theta_q_plot, plane_field, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'E_{WG}^{FF}, quartz');
hold on;
plane_field = NaN(1, length(theta_p_plot));
plane_field(1 : length(wg_field_plastic.theta)) = fliplr(Ep_total(1601, :));
plane_field(length(wg_field_plastic.theta) + 1 : end) = Ep_total(1, :);
plot(theta_s_plot, plane_field, '-.', 'LineWidth', 2.0, ...
    'DisplayName', 'E_{WG}^{FF}, plastic');
grid on;
ylim([-40 0]);
xlim([-90 90]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E^{FF}| / dB');
title('Waveguide E^{FF} @ \phi = 0 deg and 180 deg, f = 70 GHz');

figure('Position', [250 250 725 400]);
plane_field = NaN(1, length(theta_s_plot));
plane_field(1 : length(wg_field_silicon.theta)) = fliplr(Es_total(2400, :));
plane_field(length(wg_field_silicon.theta) + 1 : end) = Es_total(801, :);
plot(theta_s_plot, plane_field, 'LineWidth', 2.0, ...
    'DisplayName', 'E_{WG}^{FF}, silicon');
hold on;
plane_field = NaN(1, length(theta_q_plot));
plane_field(1 : length(wg_field_quatz.theta)) = fliplr(Eq_total(2400, :));
plane_field(length(wg_field_quatz.theta) + 1 : end) = Eq_total(801, :);
plot(theta_q_plot, plane_field, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'E_{WG}^{FF}, quartz');
hold on;
plane_field = NaN(1, length(theta_p_plot));
plane_field(1 : length(wg_field_plastic.theta)) = fliplr(Ep_total(2400, :));
plane_field(length(wg_field_plastic.theta) + 1 : end) = Ep_total(801, :);
plot(theta_s_plot, plane_field, '-.', 'LineWidth', 2.0, ...
    'DisplayName', 'E_{WG}^{FF}, plastic');
grid on;
ylim([-40 0]);
xlim([-90 90]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E^{FF}| / dB');
title('Waveguide E^{FF} @ \phi = 90 deg and 270 deg, f = 70 GHz');
