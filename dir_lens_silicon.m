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

%% PARAMETERS
% Wave parameters
wave.f = 70 * 1e9;
% Waveguide parameters
waveguide.E10 = 1;
waveguide.er = 1;
waveguide.a = 2.80 * 1e-3;
waveguide.b = 1.88 * 1e-3;
% Lens parameters
lens.er = 11.9;
theta_max = 40 * pi / 180;
% Grids parameters
Nrho = 50;
Nphi = 50;
Ntheta = 50;
% Far-field parameters
R = 1;

D_ratio = linspace(3, 15, 50);
dir_broadside = NaN(3, length(D_ratio));
for theta_max_idx = 1 : 1 : length(theta_max)
    for d_idx = 1 : 1 : length(D_ratio)
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
        lens.D = D_ratio(d_idx) * wave.wavelength;
        lens.n = sqrt(lens.er);
        lens.Z = zeta / lens.n;
        
        %% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
        waveguide_lens.TE_coef = 2 * lens.Z ./ (lens.Z + waveguide.ZTE);
        
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
        
        %% LENS PARAMETERS
        [lens, ~, lens_sph_grid, theta_inc, theta_tr] ...
            = lens_calculation(lens, cyl_grid, theta_max(theta_max_idx));
            
        %% RADIAL DISTANCE TO LENS INTERFACE AND SPHERICAL GRID INSIDE LENS
        r = lens_sph_grid(:, :, 1);
        lens_sph_grid = lens_sph_grid(:, :, 2 : 3);
        
        %% WAVE VECTOR COMPONENTS
        [k_comp, k] = wave_vector(lens.er, wave.k0, lens_sph_grid);
        
        %% WAVEGUIDE FEED
        [Ewg, Mwg] = waveguide_feed(waveguide, waveguide_lens.TE_coef, ...
            k, k_comp, r, lens_sph_grid, 'NeglectPhase');
        
        %% APERTURE CURRENT
        [J, lens, Jft] = lens_current(lens, Ewg, theta_inc, theta_tr, 1, ...
            cyl_grid, lens_sph_grid, 'FT', k0_comp);
        
        %% RADIATED FAR-FIELD
        E = farfield(wave.k0, R, sph_grid, k0_comp(:, :, 3), SGFej, Jft);
        E_total = total_field(E);
        
        %% APERTURE DIRECTIVITY, AND RADIATED POWER
        [dir, ~, rad_power] = directivity(1, E, sph_grid, R);
        dir_broadside(theta_max_idx, d_idx) = dir(1, 1);
    end
end

%% PLOT DIRECTIVITY
figure('Position', [250 250 750 400]);
for theta_max_idx = 1 : 1 : length(theta_max)
    plot(D_ratio, 10 * log10(dir_broadside(theta_max_idx, :)), ...
        'LineWidth', 2.0, 'DisplayName', ['D, \theta_{max} = ' ...
        num2str(theta_max(theta_max_idx) * 180 / pi) ' deg']);
    hold on;
end
hold off;
grid on;
xlim([min(D_ratio) max(D_ratio)]);
legend show;
legend('location', 'bestoutside');
xlabel('D / \lambda_{0}');
ylabel('D(\theta=0^{\circ},\phi=0^{\circ}) / dB');
title('Broadside Directivity @ silicon lens, f = 70 GHz');
saveas(gcf, 'figures\dir_lens_silicon.fig');

%% SAVE WORKSPACE
lens_silicon.dir_broadside = dir_broadside;
lens_silicon.D_ratio = D_ratio;
lens_silicon.theta_max = theta_max;
save('results\dir_lens_silicon.mat', 'lens_silicon');
