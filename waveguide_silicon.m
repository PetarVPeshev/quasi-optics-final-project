% Parameters to optimize:
%   Waveguide: a, b
%   Lens: D, theta_max
% Automotive radar applications:
% @ ~ 70 GHz
% 47 GHz < f < 94 GHz
% maximize BW, fc < f < 2fc, f = 1.5fc, fc = 47 GHz
% second mode: TE20
% min(a) = 3.19 mm, min(b) = 1.59 mm

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

c = physconst('LightSpeed');
zeta = 376.730313668;

Na = 101;

waveguide(Na) = struct('E10', [], 'er', [], 'a', [], 'b', [], ...
    'kx1', [], 'kz', [], 'n', [], 'ZTE', []);
waveguide_lens(Na) = struct('TE_coef', []);

%% PARAMETERS
% Wave parameters
wave.f = 70e9;
% Waveguide parameters
a = linspace(2.15, 6, Na) * 1e-3;
for a_idx = 1 : 1 : Na
    waveguide(a_idx).E10 = 1;
    waveguide(a_idx).er = 1;
    waveguide(a_idx).a = a(a_idx);
    waveguide(a_idx).b = 1.88 * 1e-3;
end
% Lens parameters
lens.er = 11.9;
% Grid parameters
Ntheta = 800;
Nphi = 3200;
% Far-field parameters
R = 1;

%% DEPENDENT PARAMETERS
% Wave parameters
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
% Waveguide parameters
for a_idx = 1 : 1 : Na
    waveguide(a_idx).kx1 = pi / waveguide(a_idx).a;
    waveguide(a_idx).kz = sqrt(wave.k0 ^ 2 - waveguide(a_idx).kx1 .^ 2);
    waveguide(a_idx).n = sqrt(waveguide(a_idx).er);
    waveguide(a_idx).ZTE = zeta * wave.k0 ./ waveguide(a_idx).kz;
end
% Lens parameters
lens.n = sqrt(lens.er);
lens.Z = zeta / lens.n;

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
for a_idx = 1 : 1 : Na
    waveguide_lens(a_idx).TE_coef = 2 * lens.Z ./ (lens.Z + waveguide(a_idx).ZTE);
end

%% COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = linspace(0, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, k] = wave_vector(lens.er, wave.k0, sph_grid);

dir_broadside = NaN(1, Na);
for a_idx = 1 : 1 : Na
    %% EQUIVALENT MAGNETIC CURRENT AND RADIATED FAR-FIELD
    % Waveguide radiated field and equivalent magnetic current
    [E, M] = waveguide_feed(waveguide(a_idx), waveguide_lens(a_idx).TE_coef, k, k_comp, R, sph_grid);
    % Waveguide total radiated magnetic field
    E_total = total_field(E);
    
    %% DIRECTIVITY
    [dir, ~, ~] = directivity(lens.er, E, sph_grid, R);
    dir_broadside(a_idx) = dir(1, 1);
end

figure('Position', [250 250 750 400]);
plot(a * 1e3, 10 * log10(dir_broadside), 'LineWidth', 2.0, ...
    'DisplayName', 'D');
hold on;
xline(2.5, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
    'DisplayName', 'min\{a\}, a = 2.5 mm');
hold on;
xline(3.5, '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0, ...
    'DisplayName', 'initial a, a = 3.5 mm');
grid on;
xlim([min(a) max(a)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('a / mm');
ylabel('D / dB');
title(['D(\theta=0^{\circ},\phi=0^{\circ}) @ f = ' ...
    num2str(wave.f * 1e-9) ' GHz, and \epsilon_{r} = ' num2str(lens.er)]);

ZTE = NaN(1, Na);
for a_idx = 1 : 1 : Na
    ZTE(a_idx) = waveguide(a_idx).ZTE;
end
figure('Position', [250 250 750 400]);
plot(a * 1e3, ZTE * 1e-3, 'LineWidth', 2.0, 'DisplayName', 'Z_{TE}');
hold on;
xline(2.5, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
    'DisplayName', 'min\{a\}, a = 2.5 mm');
hold on;
xline(3.5, '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0, ...
    'DisplayName', 'initial a, a = 3.5 mm');
grid on;
xlim([min(a) max(a)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('a / mm');
ylabel('Z_{TE} / k\Omega');
title(['Z_{TE} @ f = ' num2str(wave.f * 1e-9) ' GHz']);

P_ratio = NaN(1, Na);
for a_idx = 1 : 1 : Na
    P_ratio(a_idx) = (abs(waveguide_lens(a_idx).TE_coef) ^ 2) * waveguide(a_idx).ZTE / lens.Z;
end
figure('Position', [250 250 750 400]);
plot(a * 1e3, P_ratio, 'LineWidth', 2.0, ...
    'DisplayName', 'transmitted ratio');
hold on;
xline(2.5, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
    'DisplayName', 'min\{a\}, a = 2.5 mm');
hold on;
xline(3.5, '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0, ...
    'DisplayName', 'initial a, a = 3.5 mm');
grid on;
xlim([min(a) max(a)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('a / mm');
ylabel('P_{t}/P_{i}');
title(['TE Power Ratio Waveguide-Lens @ f = ' ...
    num2str(wave.f * 1e-9) ' GHz, and \epsilon_{r} = ' num2str(lens.er)]);

kz = NaN(1, Na);
for a_idx = 1 : 1 : Na
    kz(a_idx) = waveguide(a_idx).kz;
end
figure('Position', [250 250 750 400]);
plot(a * 1e3, kz, 'LineWidth', 2.0, 'DisplayName', 'k_{z}');
hold on;
xline(2.5, '--', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2.0, ...
    'DisplayName', 'min\{a\}, a = 2.5 mm');
hold on;
xline(3.5, '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0, ...
    'DisplayName', 'initial a, a = 3.5 mm');
grid on;
xlim([min(a) max(a)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('a / mm');
ylabel('k_{z} / rad/m');
title(['Propagation Constant @ f = ' ...
    num2str(wave.f * 1e-9) ' GHz, and \epsilon_{r} = ' num2str(lens.er)]);

%% SAVE WORKSPACE
wg_silicon = struct('a', a, 'D', dir_broadside, 'ZTE', ZTE, ...
    'P_ratio', P_ratio, 'kz', kz);
save('results\wg_silicon.mat', 'wg_silicon');
