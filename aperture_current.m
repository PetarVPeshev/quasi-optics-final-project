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
    'TM_coef', [], 'TE_coef', [], 'Efeed', [], 'Mfeed', [], 'J', []);
waveguide_lens(3) = struct('TE_coef', [], 'TE_T', []);

%% PARAMETERS
wave.f = 70e9;
waveguide.E10 = 1;
waveguide.er = 1;
waveguide.a = 3.212e-3;
waveguide.b = 1.6e-3;
lens(1).er = 11.9;
lens(2).er = 4;
lens(3).er = 2;
Nrho = 800;
Nphi = 3200;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
waveguide.kx1 = pi / waveguide.a;
waveguide.kz = sqrt(wave.k0 ^ 2 - waveguide.kx1 .^ 2);
for lens_idx = 1 : 1 : length(lens)
    lens(lens_idx).D = 6 * wave.wavelength;
    lens(lens_idx).n = sqrt(lens(lens_idx).er);
    lens(lens_idx).Z = zeta / lens(lens_idx).n;
end

%% REFRACTIVE INDECIES
waveguide.n = sqrt(waveguide.er);

%% CYLINDRICAL COORDINATE GRID INSIDE LENS
rho = linspace(eps, lens(1).D / 2, Nrho);
phi = linspace(eps, 2 * pi, Nphi);
cyl_grid = meshgrid_comb(rho, phi);

%% WAVEGUIDE TE IMPEDANCE
waveguide.ZTE = zeta * wave.k0 ./ waveguide.kz;

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
for lens_idx = 1 : 1 : length(lens)
    waveguide_lens(lens_idx).TE_coef = 2 * lens(lens_idx).Z ...
        ./ (lens(lens_idx).Z + waveguide.ZTE);
end

for lens_idx = 1 : 1 : length(lens)
    %% LENS PARAMETERS
    [lens(lens_idx), ~, sph_grid, theta_inc, theta_tr] ...
        = lens_calculation(lens(lens_idx), cyl_grid);
    
    %% RADIAL DISTANCE TO LENS INTERFACE AND SPHERICAL GRID
    r = sph_grid(:, :, 1);
    sph_grid = sph_grid(:, :, 2 : 3);

    %% WAVE VECTOR COMPONENTS
    [k_comp, k] = wave_vector(lens(lens_idx).er, wave.k0, sph_grid);

    %% WAVEGUIDE FEED
    [lens(lens_idx).Efeed, lens(lens_idx).Mfeed] ...
        = waveguide_feed(waveguide, waveguide_lens(lens_idx).TE_coef, ...
        k, k_comp, r, sph_grid, 'NeglectPhase');

    %% APERTURE CURRENT
    [J, lens(lens_idx)] = lens_current(lens(lens_idx), ...
        lens(lens_idx).Efeed, theta_inc, theta_tr, 1, cyl_grid, sph_grid);

    lens(lens_idx).theta_inc = theta_inc;
    lens(lens_idx).theta_tr = theta_tr;
    lens(lens_idx).J = J;
end

%% CARTESIAN COORDINATES
cyl_coord = zeros( [size(cyl_grid, 1, 2), 3] );
cyl_coord(:, :, 1 : 2) = cyl_grid;
cart_coord = cyl2cart_cord(cyl_coord);

%% PLOT EQUIVALENT APERTURE ELECTRIC CURRENT DENSITY
for lens_idx = 1 : 1 : length(lens)
    J_max_magn = 20 * log10( max(abs(lens(lens_idx).J), [], 'all') );
    J_norm = NaN( [size(cart_coord, 1, 2), 3] );
    J_norm(:, :, 1) = 20 * log10(abs(lens(lens_idx).J(:, :, 1))) ...
        - J_max_magn;
    J_norm(:, :, 2) = 20 * log10(abs(lens(lens_idx).J(:, :, 2))) ...
        - J_max_magn;
    J_norm(:, :, 3) = 20 * log10(abs(lens(lens_idx).J(:, :, 3))) ...
        - J_max_magn;

    figure('Position', [250 250 1050 400]);
    subplot(1, 2, 1);
    surface(cart_coord(:, :, 1) * 1e3, cart_coord(:, :, 2) * 1e3, ...
        J_norm(:, :, 1), 'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
    caxis([-40 -20]);
    view(0, 90);
    xlabel('x / mm');
    ylabel('y / mm');
    zlabel('|J_{x}| / dB');
    title('|J_{x}| / dB');
    subplot(1, 2, 2);
    surface(cart_coord(:, :, 1) * 1e3, cart_coord(:, :, 2) * 1e3, ...
        J_norm(:, :, 2), 'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
    caxis([-40 0]);
    view(0, 90);
    xlabel('x / mm');
    ylabel('y / mm');
    zlabel('|J_{y}| / dB');
    title('|J_{y}| / dB');
    sgtitle(['Untruncated Aperture J_{eq} @ f = ' num2str(wave.f * 1e-9) ...
        ' GHz, D = ' num2str(round(lens(lens_idx).D * 1e3, 2)) ' mm, ' ...
        '\epsilon_{r} = ' num2str(lens(lens_idx).er) ', a = ' ...
        num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
        num2str(round(waveguide.b * 1e3, 2)) ' mm'], 'FontSize', 17, ...
        'FontWeight', 'bold');
    saveas(gcf, ['figures\aperture_Jeq_lens_er_' ...
        num2str(lens(lens_idx).er) '.fig']);
end

%% SAVE WORKSPACE
save('results\aperture_current.mat', 'wave', 'waveguide', 'lens', ...
    'cart_coord', 'waveguide_lens', 'lens');
