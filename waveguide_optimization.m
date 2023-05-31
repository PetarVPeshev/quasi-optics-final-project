close all;
clear;
clc;

waveguide_plastic;
waveguide_quartz;
waveguide_silicon;

close all;

load('results\wg_silicon.mat');
load('results\wg_quartz.mat');
load('results\wg_plastic.mat');

figure('Position', [250 250 750 400]);
plot(wg_silicon.a * 1e3, 10 * log10(wg_silicon.D), 'LineWidth', 2.0, ...
    'DisplayName', 'D, \epsilon_{r} = 11.9');
hold on;
plot(wg_quartz.a * 1e3, 10 * log10(wg_quartz.D), 'LineWidth', 2.0, ...
    'DisplayName', 'D, \epsilon_{r} = 4');
hold on;
plot(wg_plastic.a * 1e3, 10 * log10(wg_plastic.D), 'LineWidth', 2.0, ...
    'DisplayName', 'D, \epsilon_{r} = 2');
hold on;
xline(2.5, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
    'DisplayName', 'min\{a\}, a = 2.5 mm');
hold on;
xline(3.5, '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0, ...
    'DisplayName', 'initial a, a = 3.5 mm');
grid on;
xlim([min(wg_silicon.a) max(wg_silicon.a)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('a / mm');
ylabel('D(\theta=0^{\circ},\phi=0^{\circ}) / dB');
title('Broadside Directivity @ f = 70 GHz, b = 1.88 mm');

figure('Position', [250 250 750 400]);
plot(wg_silicon.a * 1e3, wg_silicon.ZTE * 1e-3, 'LineWidth', 2.0, ...
    'DisplayName', 'Z_{TE}');
xline(2.5, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
    'DisplayName', 'min\{a\}, a = 2.5 mm');
hold on;
xline(3.5, '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0, ...
    'DisplayName', 'initial a, a = 3.5 mm');
grid on;
xlim([min(wg_silicon.a) max(wg_silicon.a)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('a / mm');
ylabel('Z_{TE} / k\Omega');
title('Waveguide Impedance @ f = 70 GHz, b = 1.88 mm');

figure('Position', [250 250 750 400]);
plot(wg_silicon.a * 1e3, wg_silicon.P_ratio, 'LineWidth', 2.0, ...
    'DisplayName', '|T^{TE}|^{2}, \epsilon_{r} = 11.9');
hold on;
plot(wg_quartz.a * 1e3, wg_quartz.P_ratio, 'LineWidth', 2.0, ...
    'DisplayName', '|T^{TE}|^{2}, \epsilon_{r} = 4');
hold on;
plot(wg_plastic.a * 1e3, wg_plastic.P_ratio, 'LineWidth', 2.0, ...
    'DisplayName', '|T^{TE}|^{2}, \epsilon_{r} = 2');
hold on;
xline(2.5, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
    'DisplayName', 'min\{a\}, a = 2.5 mm');
hold on;
xline(3.5, '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0, ...
    'DisplayName', 'initial a, a = 3.5 mm');
grid on;
xlim([min(wg_silicon.a) max(wg_silicon.a)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('a / mm');
ylabel('P_{t}/P_{i}');
title('TE Transmitted Power Ratio Waveguide-Lens @ f = 70 GHz, b = 1.88 mm');

% figure('Position', [250 250 750 400]);
% plot(wg_silicon.a * 1e3, wg_silicon.kz, 'LineWidth', 2.0, ...
%     'DisplayName', 'k_{z}');
% xline(3.19, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
%     'DisplayName', 'min\{a\}');
% hold on;
% xline(3.5, '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0, ...
%     'DisplayName', 'initial a, a = 3.5 mm');
% grid on;
% xlim([min(wg_silicon.a) max(wg_silicon.a)] * 1e3);
% legend show;
% legend('location', 'bestoutside');
% xlabel('a / mm');
% ylabel('k_{z} / rad/m');
% title('Propagation Constant @ f = 70 GHz, b = 1.88 mm');

