% Parameters to optimize:
%   Waveguide: a, b
%   Lens: D, theta_max
% Automotive radar applications:
% REC. ITU-R M.2057-1
% @ ~ 70 GHz, Type A radar
% BW = 1 GHz
% Allocate BW in waveguide BW = 20 GHz
% 65 GHz < f < 75 GHz
% second mode: TE01
% min(a) = 2.5 mm, max(b) = 2 mm

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
waveguide(7) = struct('E10', [], 'er', [], 'a', [], 'b', [], 'kx1', [], ...
    'kz', [], 'n', [], 'ZTE', []);
lens = struct('er', [], 'n', [], 'Z', [], 'TE_coef', [], 'k', [], ...
    'k_comp', [], 'M', [], 'E', [], 'E_total', [], 'dir', [], ...
    'dir_broadside', []);
waveguide_field(7) = struct('lens_silicon', [], 'lens_quartz', [], ...
    'lens_plastic', []);

%% PARAMETERS
% Wave parameters
wave.f = 70e9;
% Waveguide parameters
for a_idx = 1 : 1 : length(waveguide)
    waveguide(a_idx).E10 = 1;
    waveguide(a_idx).er = 1;
    waveguide(a_idx).a = (20 + 0.2 * a_idx) * 1e-3;
    waveguide(a_idx).b = 2e-3;
end
% Lens parameters
for wg_idx = 1 : 1 : length(waveguide)
    waveguide_field(wg_idx).lens_silicon = lens;
    waveguide_field(wg_idx).lens_silicon.er = 11.9;
    waveguide_field(wg_idx).lens_quartz = lens;
    waveguide_field(wg_idx).lens_quartz.er = 4;
    waveguide_field(wg_idx).lens_plastic = lens;
    waveguide_field(wg_idx).lens_plastic.er = 2;
end
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
for wg_idx = 1 : 1 : length(waveguide)
    waveguide(wg_idx).kx1 = pi ./ waveguide(wg_idx).a;
    waveguide(wg_idx).kz = sqrt(wave.k0 ^ 2 - waveguide(wg_idx).kx1 .^ 2);
    waveguide(wg_idx).n = sqrt(waveguide(wg_idx).er);
    waveguide(wg_idx).ZTE = zeta * wave.k0 ./ waveguide(wg_idx).kz;
end
% Lens parameters
for wg_idx = 1 : 1 : length(waveguide)
    waveguide_field(wg_idx).lens_silicon.n = sqrt(waveguide_field(wg_idx).lens_silicon.er);
    waveguide_field(wg_idx).lens_silicon.Z = zeta / waveguide_field(wg_idx).lens_silicon.n;
    waveguide_field(wg_idx).lens_quartz.n = sqrt(waveguide_field(wg_idx).lens_quartz.er);
    waveguide_field(wg_idx).lens_quartz.Z = zeta / waveguide_field(wg_idx).lens_quartz.n;
    waveguide_field(wg_idx).lens_plastic.n = sqrt(waveguide_field(wg_idx).lens_plastic.er);
    waveguide_field(wg_idx).lens_plastic.Z = zeta / waveguide_field(wg_idx).lens_plastic.n;
end

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
for wg_idx = 1 : 1 : length(waveguide)
    waveguide_field(wg_idx).lens_silicon.TE_coef = 2 * waveguide_field(wg_idx).lens_silicon.Z ...
        ./ (waveguide_field(wg_idx).lens_silicon.Z + waveguide(wg_idx).ZTE);
    waveguide_field(wg_idx).lens_quartz.TE_coef = 2 * waveguide_field(wg_idx).lens_quartz.Z ...
        ./ (waveguide_field(wg_idx).lens_quartz.Z + waveguide(wg_idx).ZTE);
    waveguide_field(wg_idx).lens_plastic.TE_coef = 2 * waveguide_field(wg_idx).lens_plastic.Z ...
        ./ (waveguide_field(wg_idx).lens_plastic.Z + waveguide(wg_idx).ZTE);
end

%% COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = linspace(0, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

for wg_idx = 1 : 1 : length(waveguide)
    %% WAVE VECTOR COMPONENTS
    [waveguide_field(wg_idx).lens_silicon.k_comp, waveguide_field(wg_idx).lens_silicon.k] ...
        = wave_vector(waveguide_field(wg_idx).lens_silicon.er, wave.k0, sph_grid);
    [waveguide_field(wg_idx).lens_quartz.k_comp, waveguide_field(wg_idx).lens_quartz.k] ...
        = wave_vector(waveguide_field(wg_idx).lens_quartz.er, wave.k0, sph_grid);
    [waveguide_field(wg_idx).lens_plastic.k_comp, waveguide_field(wg_idx).lens_plastic.k] ...
        = wave_vector(waveguide_field(wg_idx).lens_plastic.er, wave.k0, sph_grid);

    %% EQUIVALENT MAGNETIC CURRENT AND RADIATED FAR-FIELD
    % Waveguide radiated field and equivalent magnetic current
    [waveguide_field(wg_idx).lens_silicon.E, waveguide_field(wg_idx).lens_silicon.M] ...
        = waveguide_feed(waveguide, waveguide_field(wg_idx).lens_silicon.TE_coef, ...
        waveguide_field(wg_idx).lens_silicon.k, waveguide_field(wg_idx).lens_silicon.k_comp, ...
        R, sph_grid);
    [waveguide_field(wg_idx).lens_quartz.E, waveguide_field(wg_idx).lens_quartz.M] ...
        = waveguide_feed(waveguide, waveguide_field(wg_idx).lens_quartz.TE_coef, ...
        waveguide_field(wg_idx).lens_quartz.k, waveguide_field(wg_idx).lens_quartz.k_comp, ...
        R, sph_grid);
    [waveguide_field(wg_idx).lens_plastic.E, waveguide_field(wg_idx).lens_plastic.M] ...
        = waveguide_feed(waveguide, waveguide_field(wg_idx).lens_plastic.TE_coef, ...
        waveguide_field(wg_idx).lens_plastic.k, waveguide_field(wg_idx).lens_plastic.k_comp, ...
        R, sph_grid);
    % Waveguide total radiated magnetic field
    waveguide_field(wg_idx).lens_silicon.E_total = total_field(waveguide_field(wg_idx).lens_silicon.E);
    waveguide_field(wg_idx).lens_quartz.E_total = total_field(waveguide_field(wg_idx).lens_quartz.E);
    waveguide_field(wg_idx).lens_plastic.E_total = total_field(waveguide_field(wg_idx).lens_plastic.E);

    %% DIRECTIVITY
    waveguide_field(wg_idx).lens_silicon.dir = directivity(waveguide_field(wg_idx).lens_silicon.er, ...
        waveguide_field(wg_idx).lens_silicon.E, sph_grid, R);
    waveguide_field(wg_idx).lens_silicon.dir_broadside = waveguide_field(wg_idx).lens_silicon.dir(1, 1);
    waveguide_field(wg_idx).lens_quartz.dir = directivity(waveguide_field(wg_idx).lens_quartz.er, ...
        waveguide_field(wg_idx).lens_quartz.E, sph_grid, R);
    waveguide_field(wg_idx).lens_quartz.dir_broadside = waveguide_field(wg_idx).lens_quartz.dir(1, 1);
    waveguide_field(wg_idx).lens_plastic.dir = directivity(waveguide_field(wg_idx).lens_plastic.er, ...
        waveguide_field(wg_idx).lens_plastic.E, sph_grid, R);
    waveguide_field(wg_idx).lens_plastic.dir_broadside = waveguide_field(wg_idx).lens_plastic.dir(1, 1);
end

%% PLOT DIRECTIVITY 
