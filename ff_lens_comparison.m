close all;
clear;
clc;

ff_lens_silicon;
ff_lens_quartz;
ff_lens_plastic;

close all;

load('results\lens_silicon.mat');
load('results\lens_quartz.mat');
load('results\lens_plastic.mat');

Es_total = norm_magnitude(total_field(lens_silicon.E), 'dB');
Ep_total = norm_magnitude(total_field(lens_plastic.E), 'dB');
Eq_total = norm_magnitude(total_field(lens_quartz.E), 'dB');

theta_s_plot = NaN(1, 2 * length(lens_silicon.theta));
theta_s_plot(1 : length(lens_silicon.theta)) = - fliplr(lens_silicon.theta) * 180 / pi;
theta_s_plot(length(lens_silicon.theta) + 1 : end) = lens_silicon.theta * 180 / pi;
theta_q_plot = NaN(1, 2 * length(lens_quartz.theta));
theta_q_plot(1 : length(lens_quartz.theta)) = - fliplr(lens_quartz.theta) * 180 / pi;
theta_q_plot(length(lens_quartz.theta) + 1 : end) = lens_quartz.theta * 180 / pi;
theta_p_plot = NaN(1, 2 * length(lens_plastic.theta));
theta_p_plot(1 : length(lens_plastic.theta)) = - fliplr(lens_plastic.theta) * 180 / pi;
theta_p_plot(length(lens_plastic.theta) + 1 : end) = lens_plastic.theta * 180 / pi;

figure('Position', [250 250 725 400]);
plane_field = NaN(1, length(theta_s_plot));
plane_field(1 : length(lens_silicon.theta)) = fliplr(Es_total(51, :));
plane_field(length(lens_silicon.theta) + 1 : end) = Es_total(1, :);
plot(theta_s_plot, plane_field, 'LineWidth', 2.0, ...
    'DisplayName', 'E^{FF}, silicon');
hold on;
plane_field = NaN(1, length(theta_q_plot));
plane_field(1 : length(lens_quartz.theta)) = fliplr(Eq_total(51, :));
plane_field(length(lens_quartz.theta) + 1 : end) = Eq_total(1, :);
plot(theta_q_plot, plane_field, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'E^{FF}, quartz');
hold on;
plane_field = NaN(1, length(theta_p_plot));
plane_field(1 : length(lens_plastic.theta)) = fliplr(Ep_total(51, :));
plane_field(length(lens_plastic.theta) + 1 : end) = Ep_total(1, :);
plot(theta_s_plot, plane_field, '-.', 'LineWidth', 2.0, ...
    'DisplayName', 'E^{FF}, plastic');
grid on;
ylim([-35 0]);
xlim([-30 30]);
xticks(-30 : 5 : 30);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E^{FF}| / dB');
title('Lens E^{FF} @ \phi = 0 deg and 180 deg, f = 70 GHz');

figure('Position', [250 250 725 400]);
plane_field = NaN(1, length(theta_s_plot));
plane_field(1 : length(lens_silicon.theta)) = fliplr(Es_total(75, :));
plane_field(length(lens_silicon.theta) + 1 : end) = Es_total(26, :);
plot(theta_s_plot, plane_field, 'LineWidth', 2.0, ...
    'DisplayName', 'E^{FF}, silicon');
hold on;
plane_field = NaN(1, length(theta_q_plot));
plane_field(1 : length(lens_quartz.theta)) = fliplr(Eq_total(75, :));
plane_field(length(lens_quartz.theta) + 1 : end) = Eq_total(26, :);
plot(theta_q_plot, plane_field, '--', 'LineWidth', 2.0, ...
    'DisplayName', 'E^{FF}, quartz');
hold on;
plane_field = NaN(1, length(theta_p_plot));
plane_field(1 : length(lens_plastic.theta)) = fliplr(Ep_total(75, :));
plane_field(length(lens_plastic.theta) + 1 : end) = Ep_total(26, :);
plot(theta_s_plot, plane_field, '-.', 'LineWidth', 2.0, ...
    'DisplayName', 'E^{FF}, plastic');
grid on;
ylim([-35 0]);
xlim([-30 30]);
xticks(-30 : 5 : 30);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E^{FF}| / dB');
title('Lens E^{FF} @ \phi = 90 deg and 270 deg, f = 70 GHz');
