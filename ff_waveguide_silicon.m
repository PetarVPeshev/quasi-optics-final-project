% Parameters to optimize:
%   Waveguide: a, b
%   Lens: D, theta_max
% Automotive radar applications:
% REC. ITU-R M.2057-1
% @ ~ 70 GHz, Type A radar
% BW = 1 GHz
% Allocate BW in waveguide BW = 20 GHz
% 60 GHz < f < 80 GHz
% second mode: TE01
% min(a) = 2.5 mm, max(b) = 1.87 mm

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

Na = 2;
Nb = 2;

waveguide(Na, Nb) = struct('E10', [], 'er', [], 'a', [], 'b', [], ...
    'kx1', [], 'kz', [], 'n', [], 'ZTE', []);
waveguide_lens(Na, Nb) = struct('TE_coef', []);

%% PARAMETERS
% Wave parameters
wave.f = 70e9;
% Waveguide parameters
a = linspace(2.5, 5, Na) * 1e-3;
b = linspace(1.25, 2.5, Nb) * 1e-3;
for a_idx = 1 : 1 : Na
    for b_idx = 1 : 1 : Nb
        waveguide(a_idx, b_idx).E10 = 1;
        waveguide(a_idx, b_idx).er = 1;
        waveguide(a_idx, b_idx).a = a(a_idx);
        waveguide(a_idx, b_idx).b = b(b_idx);
    end
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
    for b_idx = 1 : 1 : Nb
        waveguide(a_idx, b_idx).kx1 = pi / waveguide(a_idx, b_idx).a;
        waveguide(a_idx, b_idx).kz = sqrt(wave.k0 ^ 2 - waveguide(a_idx, b_idx).kx1 .^ 2);
        waveguide(a_idx, b_idx).n = sqrt(waveguide(a_idx, b_idx).er);
        waveguide(a_idx, b_idx).ZTE = zeta * wave.k0 ./ waveguide(a_idx, b_idx).kz;
    end
end
% Lens parameters
lens.n = sqrt(lens.er);
lens.Z = zeta / lens.n;

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
for a_idx = 1 : 1 : Na
    for b_idx = 1 : 1 : Nb
        waveguide_lens(a_idx, b_idx).TE_coef = 2 * lens.Z ...
            ./ (lens.Z + waveguide(a_idx, b_idx).ZTE);
    end
end

%% COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = linspace(0, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, k] = wave_vector(lens.er, wave.k0, sph_grid);

figure('Position', [250 250 750 400]);
dir_broadside = NaN(Na, Nb);
for a_idx = 1 : 1 : Na
    for b_idx = 1 : 1 : Nb
        %% EQUIVALENT MAGNETIC CURRENT AND RADIATED FAR-FIELD
        % Waveguide radiated field and equivalent magnetic current
        [E, M] = waveguide_feed(waveguide(a_idx, b_idx), ...
            waveguide_lens(a_idx, b_idx).TE_coef, k, k_comp, R, sph_grid);
        % Waveguide total radiated magnetic field
        E_total = total_field(E);
        
        %% DIRECTIVITY
        dir = directivity(lens.er, E, sph_grid, R);
        dir_broadside(a_idx, b_idx) = dir(1, 1);
        
        E_plot = NaN(1, 2 * length(theta));
        E_plot(1 : length(theta)) = fliplr(E_total(1601, :));
        E_plot(length(theta) + 1 : end) = E_total(1, :);
        E_plot = 20 * log10(E_plot) - 20 * log10(max(E_plot, [], 'all'));
        theta_plot = NaN(1, 2 * length(theta));
        theta_plot(1 : length(theta)) = - fliplr(theta * 180 / pi);
        theta_plot(length(theta) + 1 : end) = theta * 180 / pi;
        plot(theta_plot, E_plot, 'LineWidth', 2.0, ...
            'DisplayName', ['a = ' num2str(waveguide(a_idx, b_idx).a * 1e3) ...
            ' mm, b = ' num2str(waveguide(a_idx, b_idx).b * 1e3) ' mm']);
        hold on;
    end
end
hold off;
grid on;
xlim([-90 90]);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E^{FF}| / dB');
title(['|E_{WG}^{FF}| @ f = ' num2str(wave.f * 1e-9) ' GHz, and ' ...
    '\epsilon_{r} = ' num2str(lens.er)]);

[A, B] = meshgrid(a, b);
figure();
surface(A, B, 10 * log10(dir_broadside), 'LineStyle', 'none');
grid on;
colormap('jet');
colorbar;
