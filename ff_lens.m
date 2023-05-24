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
lens(3) = struct('er', [], 'D', [], 'n', [], 'Z', [], 'R', [], 'e', [], ...
    'theta_crit', [], 'theta_max', [], 'r_min', [], 'a', [], 'c', [], ...
    'b', [], 'theta_inc', [], 'theta_tr', [], ...
    'TM_coef', [], 'TE_coef', [], 'Efeed', [], 'Mfeed', [], ...
    'J', [], 'Jft', [], 'E', [], 'E_total', [], ...
    'dir_db', [], 'dir_broadside_db', []);
waveguide_lens(3) = struct('TE_coef', [], 'TE_T', []);

%% PARAMETERS
% Wave parameters
wave.f = 70e9;
% Waveguide parameters
waveguide.E10 = 1;
waveguide.er = 1;
waveguide.a = 3.212e-3;
waveguide.b = 1.6e-3;
% Lens parameters
lens(1).er = 11.9;
lens(2).er = 4;
lens(3).er = 2;
% Grids parameters
Nrho = 200;
Nphi = 200;
Ntheta = 200;
% Far-field parameters
R = 1;

%% DEPENDENT PARAMETERS
% Wave parameters
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
% Waveguide parameters
waveguide.kx1 = pi / waveguide.a;
waveguide.kz = sqrt(wave.k0 ^ 2 - waveguide.kx1 .^ 2);
waveguide.ZTE = zeta * wave.k0 ./ waveguide.kz;
waveguide.n = sqrt(waveguide.er);
% Lens parameters
for lens_idx = 1 : 1 : length(lens)
    lens(lens_idx).D = 6 * wave.wavelength;
    lens(lens_idx).n = sqrt(lens(lens_idx).er);
    lens(lens_idx).Z = zeta / lens(lens_idx).n;
end

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
for lens_idx = 1 : 1 : length(lens)
    waveguide_lens(lens_idx).TE_coef = 2 * lens(lens_idx).Z ...
        ./ (lens(lens_idx).Z + waveguide.ZTE);
end

%% CYLINDRICAL COORDINATE GRID INSIDE LENS
rho = linspace(eps, lens(1).D / 2, Nrho);
phi = linspace(0, 2 * pi, Nphi);
cyl_grid = meshgrid_comb(rho, phi);

%% SPHERICAL COORDINATE GRID OUTSIDE LENS
theta = linspace(0,  pi / 2 - 0.1 * pi / 180, Ntheta);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k0_comp, ~] = wave_vector(1, wave.k0, sph_grid);

%% SPECTRAL GREEN'S FUNCTION
SGFej = dyadic_sgf(1, wave.k0, k0_comp, 'E', 'J');

for lens_idx = 1 : 1 : length(lens)
    %% LENS PARAMETERS
    [lens(lens_idx), ~, lens_sph_grid, theta_inc, theta_tr] ...
        = lens_calculation(lens(lens_idx), cyl_grid);
    
    %% RADIAL DISTANCE TO LENS INTERFACE AND SPHERICAL GRID INSIDE LENS
    r = lens_sph_grid(:, :, 1);
    lens_sph_grid = lens_sph_grid(:, :, 2 : 3);

    %% WAVE VECTOR COMPONENTS
    [k_comp, k] = wave_vector(lens(lens_idx).er, wave.k0, lens_sph_grid);

    %% WAVEGUIDE FEED
    [lens(lens_idx).Efeed, lens(lens_idx).Mfeed] ...
        = waveguide_feed(waveguide, waveguide_lens(lens_idx).TE_coef, ...
        k, k_comp, r, lens_sph_grid, 'NeglectPhase');

    %% APERTURE CURRENT
    [J, lens(lens_idx), Jft] ...
        = lens_current(lens(lens_idx), lens(lens_idx).Efeed, theta_inc, ...
        theta_tr, 1, cyl_grid, lens_sph_grid, 'FT', k0_comp);

    %% RADIATED FAR-FIELD
    lens(lens_idx).E ...
        = farfield(wave.k0, R, sph_grid, k0_comp(:, :, 3), SGFej, Jft);
    lens(lens_idx).E_total = total_field(lens(lens_idx).E);

    %% APERTURE DIRECTIVITY, AND RADIATED POWER
    [dir, ~, lens(lens_idx).rad_power] ...
        = directivity(1, lens(lens_idx).E, sph_grid, R);
    lens(lens_idx).dir_db = 10 * log10(dir);
    lens(lens_idx).dir_broadside_db = lens(lens_idx).dir_db(1, 1);

    lens(lens_idx).theta_inc = theta_inc;
    lens(lens_idx).theta_tr = theta_tr;
    lens(lens_idx).J = J;
    lens(lens_idx).Jft = Jft;
end

%% UV COORDINATES
uv_grid = uv_repr(sph_grid);

for lens_idx = 1 : 1 : length(lens)
    figure('Position', [250 250 1050 400]);
    subplot(1, 2, 1);
    surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
        norm_magnitude(lens(lens_idx).E_total, 'dB'), 'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
    xlim([-1 1]);
    ylim([-1 1]);
    caxis([-40 0]);
    view(0, 90);
    xticks(-1 : 1 : 1);
    yticks(-1 : 0.5 : 1);
    xlabel('U');
    ylabel('V');
    zlabel('|E| / dB');
    subplot(1, 2, 2);
    surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
        norm_magnitude(lens(lens_idx).E_total, 'dB'), 'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
    xlim([-1 1]);
    ylim([-1 1]);
    caxis([-40 0]);
    view(-37.5, 30);
    xticks(-1 : 1 : 1);
    yticks(-1 : 0.5 : 1);
    xlabel('U');
    ylabel('V');
    zlabel('|E| / dB');
    sgtitle(['E^{FF} @ f = ' num2str(wave.f * 1e-9) ...
        ' GHz, D = ' num2str(round(lens(lens_idx).D * 1e3, 2)) ' mm, ' ...
        '\epsilon_{r} = ' num2str(lens(lens_idx).er) ', a = ' ...
        num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
        num2str(round(waveguide.b * 1e3, 2)) ' mm'], 'FontSize', 17, ...
        'FontWeight', 'bold');
    saveas(gcf, ['figures\e_ff_lens_er_' num2str(lens(lens_idx).er) ...
        '.fig']);
end

%% SAVE WORKSPACE
save('results\ff_lens.mat', 'wave', 'waveguide', 'lens', ...
    'waveguide_lens', 'sph_grid', 'R');
