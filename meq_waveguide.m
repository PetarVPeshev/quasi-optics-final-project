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

wave = struct('f', [], 'wavelength', [], 'k0', []);
waveguide = struct('E10', [], 'er', [], 'a', [], 'b', [], 'kx1', [], ...
    'kz', [], 'n', [], 'ZTE', []);
lens(3) = struct('er', [], 'n', [], 'Z', []);
waveguide_lens(3) = struct('TE_coef', []);
waveguide_meq(3) = struct('k', [], 'k_comp', [], 'M', []);

%% PARAMETERS
% Wave parameters
wave.f = 70e9;
% Waveguide parameters
waveguide.E10 = 1;
waveguide.er = 1;
waveguide.a = 2.50e-3;  % 3.212
waveguide.b = 1.87e-3;  % 1.6
% Lens parameters
lens(1).er = 11.9;
lens(2).er = 4;
lens(3).er = 2;
% Grid parameters
Ntheta = 800;
% Far-field parameters
R = 1;

%% DEPENDENT PARAMETERS
% Wave parameters
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
% Waveguide parameters
waveguide.kx1 = pi / waveguide.a;
waveguide.kz = sqrt(wave.k0 ^ 2 - waveguide.kx1 .^ 2);
waveguide.n = sqrt(waveguide.er);
waveguide.ZTE = zeta * wave.k0 ./ waveguide.kz;
% Lens parameters
for lens_idx = 1 : 1 : length(lens)
    lens(lens_idx).n = sqrt(lens(lens_idx).er);
    lens(lens_idx).Z = zeta / lens(lens_idx).n;
end

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
for lens_idx = 1 : 1 : length(lens)
    waveguide_lens(lens_idx).TE_coef = 2 * lens(lens_idx).Z ...
        ./ (lens(lens_idx).Z + waveguide.ZTE);
end

%% COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = [0 180] * pi / 180;
sph_grid = meshgrid_comb(theta, phi);

for lens_idx = 1 : 1 : length(lens)
    %% WAVE VECTOR COMPONENTS
    [waveguide_meq(lens_idx).k_comp, waveguide_meq(lens_idx).k] ...
        = wave_vector(lens(lens_idx).er, wave.k0, sph_grid);

    %% EQUIVALENT MAGNETIC CURRENT AND RADIATED FAR-FIELD
    % Waveguide radiated field and equivalent magnetic current
    [~, waveguide_meq(lens_idx).M] = waveguide_feed(waveguide, waveguide_lens(lens_idx).TE_coef, ...
        waveguide_meq(lens_idx).k, waveguide_meq(lens_idx).k_comp, R, sph_grid);
end

%% PLOT MAGNETIC CURRENT DISTRIBUTION @ E-PLANE
theta_plot = NaN(1, length(theta) * 2);
theta_plot(1 : length(theta)) = - fliplr(theta) * 180 / pi;
theta_plot(length(theta) + 1 : end) = theta * 180 / pi;
for lens_idx = 1 : 1 : length(lens)
    plane_field = NaN(1, length(theta) * 2);
    plane_field(1 : length(theta)) = fliplr(waveguide_meq(lens_idx).M(2, :, 1));
    plane_field(length(theta) + 1 : end) = waveguide_meq(lens_idx).M(1, :, 1);
    figure('Position', [250 250 750 400]);
    plot(theta_plot, abs(plane_field), 'LineWidth', 2.0, ...
        'DisplayName', ['M_{eq}^{x}, a = ' num2str(waveguide.a * 1e3) ...
        ' mm, b = ' num2str(waveguide.b * 1e3) ' mm']);
    grid on;
    xlim([-90 90]);
    legend show;
    legend('location', 'bestoutside');
    xlabel('\theta / deg');
    ylabel('|M_{eq}|');
    title(['M_{eq}^{x} @ f = ' num2str(wave.f * 1e-9) ' GHz, and ' ...
        '\epsilon_{r} = ' num2str(lens(lens_idx).er)]);
end
