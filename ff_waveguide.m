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
waveguide_field(3) = struct('k', [], 'k_comp', [], 'M', [], 'E', [], ...
    'E_total', [], 'dir', [], 'dir_broadside', []);

%% PARAMETERS
wave.f = 70e9;
waveguide.E10 = 1;
waveguide.er = 1;
waveguide.a = 3.212e-3;
waveguide.b = 1.6e-3;
lens(1).er = 11.9;
lens(2).er = 4;
lens(3).er = 2;
Ntheta = 800;
Nphi = 3200;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
waveguide.kx1 = pi / waveguide.a;
waveguide.kz = sqrt(wave.k0 ^ 2 - waveguide.kx1 .^ 2);
for lens_idx = 1 : 1 : length(lens)
    lens(lens_idx).n = sqrt(lens(lens_idx).er);
    lens(lens_idx).Z = zeta / lens(lens_idx).n;
end

%% COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = linspace(0, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% REFRACTIVE INDECIES
waveguide.n = sqrt(waveguide.er);

%% WAVEGUIDE TE IMPEDANCE
waveguide.ZTE = zeta * wave.k0 ./ waveguide.kz;

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
for lens_idx = 1 : 1 : length(lens)
    waveguide_lens(lens_idx).TE_coef = 2 * lens(lens_idx).Z ...
        ./ (lens(lens_idx).Z + waveguide.ZTE);
end

for lens_idx = 1 : 1 : length(lens)
    %% WAVE VECTOR COMPONENTS
    [waveguide_field(lens_idx).k_comp, waveguide_field(lens_idx).k] ...
        = wave_vector(lens(lens_idx).er, wave.k0, sph_grid);

    %% EQUIVALENT MAGNETIC CURRENT AND RADIATED FAR-FIELD
    % Waveguide radiated field and equivalent magnetic current
    [waveguide_field(lens_idx).E, waveguide_field(lens_idx).M] ...
        = waveguide_feed(waveguide, waveguide_lens(lens_idx).TE_coef, ...
        waveguide_field(lens_idx).k, waveguide_field(lens_idx).k_comp, ...
        R, sph_grid);
    % Waveguide total radiated magnetic field
    waveguide_field(lens_idx).E_total ...
        = total_field(waveguide_field(lens_idx).E);

    %% DIRECTIVITY
    waveguide_field(lens_idx).dir = directivity(lens(lens_idx).er, ...
        waveguide_field(lens_idx).E, sph_grid, R);
    waveguide_field(lens_idx).dir_broadside ...
        = waveguide_field(lens_idx).dir(1, 1);
end

%% UV COORDINATES
uv_grid = uv_repr(sph_grid);

for lens_idx = 1 : 1 : length(lens)
    %% PLOT MAGNETIC CURRENT DENSITY
    figure('Position', [250 250 1050 400]);
    subplot(1, 2, 1);
    surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
        norm_magnitude(waveguide_field(lens_idx).M(:, :, 1), 'dB'), ...
        'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
    xlim([-1 1]);
    ylim([-1 1]);
    caxis([-10 0]);
    view(0, 90);
    xticks(-1 : 1 : 1);
    yticks(-1 : 0.5 : 1);
    xlabel('U');
    ylabel('V');
    zlabel('|M_{eq}| / dB');
    subplot(1, 2, 2);
    surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
        norm_magnitude(waveguide_field(lens_idx).M(:, :, 1), 'dB'), ...
        'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
    xlim([-1 1]);
    ylim([-1 1]);
    caxis([-10 0]);
    view(-37.5, 30);
    xticks(-1 : 1 : 1);
    yticks(-1 : 0.5 : 1);
    xlabel('U');
    ylabel('V');
    zlabel('|M_{eq}| / dB');
    sgtitle(['Waveguide M_{eq} @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
        '\epsilon_{r,lens} = ' num2str(lens(lens_idx).er) ', a = ' ...
        num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
        num2str(round(waveguide.b * 1e3, 2)) ' mm'], 'FontSize', ...
        17, 'FontWeight', 'bold');
    saveas(gcf, ['figures\waveguide_Meq_lens_er_' ...
        num2str(lens(lens_idx).er) '.fig']);

    %% PLOT APERTURE ELECTRIC FAR-FIELD
    figure('Position', [250 250 1050 400]);
    subplot(1, 2, 1);
    surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
        norm_magnitude(waveguide_field(lens_idx).E_total, 'dB'), ...
        'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
    xlim([-1 1]);
    ylim([-1 1]);
    caxis([-20 0]);
    view(0, 90);
    xticks(-1 : 1 : 1);
    yticks(-1 : 0.5 : 1);
    xlabel('U');
    ylabel('V');
    zlabel('|E| / dB');
    subplot(1, 2, 2);
    surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
        norm_magnitude(waveguide_field(lens_idx).E_total, 'dB'), ...
        'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
    xlim([-1 1]);
    ylim([-1 1]);
    caxis([-20 0]);
    view(-37.5, 30);
    xticks(-1 : 1 : 1);
    yticks(-1 : 0.5 : 1);
    xlabel('U');
    ylabel('V');
    zlabel('|E| / dB');
    sgtitle(['Waveguide E^{FF} @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
        '\epsilon_{r,lens} = ' num2str(lens(lens_idx).er) ', a = ' ...
        num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
        num2str(round(waveguide.b * 1e3, 2)) ' mm'], 'FontSize', ...
        17, 'FontWeight', 'bold');
    saveas(gcf, ['figures\waveguide_pattern_lens_er_' ...
        num2str(lens(lens_idx).er) '.fig']);
end

%% SAVE WORKSPACE
save('results\ff_waveguide.mat', 'wave', 'R', 'waveguide', 'lens', ...
    'sph_grid', 'waveguide_lens', 'waveguide_field');
