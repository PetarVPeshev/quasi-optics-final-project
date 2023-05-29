close all;
clear;
clc;

ff_waveguide_plastic;
ff_waveguide_quartz;
ff_waveguide_silicon;

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
xline(3.19, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
    'DisplayName', 'min\{a\}');
grid on;
xlim([min(wg_silicon.a) max(wg_silicon.a)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('a / mm');
ylabel('D(\theta=0^{\circ},\phi=0^{\circ}) / dB');
title('Broadside Directivity @ f = 70 GHz');

figure('Position', [250 250 750 400]);
plot(wg_silicon.a * 1e3, wg_silicon.ZTE * 1e-3, 'LineWidth', 2.0, ...
    'DisplayName', 'Z_{TE}');
xline(3.19, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
    'DisplayName', 'min\{a\}');
grid on;
xlim([min(wg_silicon.a) max(wg_silicon.a)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('a / mm');
ylabel('Z_{TE} / k\Omega');
title('Waveguide Impedance @ f = 70 GHz');

figure('Position', [250 250 750 400]);
plot(wg_silicon.a * 1e3, wg_silicon.P_ratio, 'LineWidth', 2.0, ...
    'DisplayName', 'power ratio, \epsilon_{r} = 11.9');
hold on;
plot(wg_quartz.a * 1e3, wg_quartz.P_ratio, 'LineWidth', 2.0, ...
    'DisplayName', 'power ratio, \epsilon_{r} = 4');
hold on;
plot(wg_plastic.a * 1e3, wg_plastic.P_ratio, 'LineWidth', 2.0, ...
    'DisplayName', 'power ratio, \epsilon_{r} = 2');
hold on;
xline(3.19, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
    'DisplayName', 'min\{a\}');
grid on;
xlim([min(wg_silicon.a) max(wg_silicon.a)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('a / mm');
ylabel('P_{t}/P_{i}');
title('TE Power Ratio Waveguide-Lens @ f = 70 GHz');
