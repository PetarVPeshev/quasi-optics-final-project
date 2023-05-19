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

%% PARAMETERS
wave.f = 70e9;
waveguide.E10 = 1;
waveguide.er = 1;
waveguide.a = 3.212e-3;
waveguide.b = 1.6e-3;
lens.er = [11.9 4 2];
Ntheta = 1200;
Nphi = 4800;
R = 1;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;

%% COORDINATE GRID
theta = linspace(eps, pi / 2, Ntheta);
phi = linspace(0, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% REFRACTIVE INDECIES
waveguide.n = sqrt(waveguide.er);
lens.n = sqrt(lens.er);

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
T = 2 * waveguide.n ./ (waveguide.n + lens.n);

% %% WAVEGUIDE PROPAGATION CONSTANT
% [waveguide.k_comp, ~] = wave_vector(waveguide.er, wave.k0, sph_grid);

k = NaN(1, length(lens.er));
k_comp = NaN( [size(sph_grid, 1, 2), 3, length(lens.er)] );
M = zeros( [size(sph_grid, 1, 2), 3, length(lens.er)] );
SGFem = NaN( [size(sph_grid, 1, 2), 3, 3, length(lens.er)] );
Eff = NaN( [size(sph_grid, 1, 2), 3, length(lens.er)] );
Eff_total = NaN( [size(sph_grid, 1, 2), length(lens.er)] );
for media_idx = 1 : 1 : length(lens.er)
    %% WAVE VECTOR COMPONENTS
    [k_comp(:, :, :, media_idx), k(media_idx)] ...
        = wave_vector(lens.er(media_idx), wave.k0, sph_grid);

    %% FOURIER TRANSFORM OF EQUIVALENT MAGNETIC CURRENT
    M(:, :, 1, media_idx) = 4 * pi * T(media_idx) * waveguide.E10 ...
        * waveguide.b * cos(waveguide.a * k_comp(:, :, 1, media_idx) / 2) ...
        .* sinc(waveguide.b * k_comp(:, :, 2, media_idx) / 2) ...
        ./ (waveguide.a * (k_comp(:, :, 1, media_idx) .^ 2 ...
        - (pi / waveguide.a) ^ 2));
%     M(:, :, 1, media_idx) = 4 * pi * T(media_idx) * waveguide.E10 ...
%         * waveguide.b * cos(waveguide.a * waveguide.k_comp(:, :, 1) / 2) ...
%         .* sinc(waveguide.b * waveguide.k_comp(:, :, 2) / 2) ...
%         ./ (waveguide.a * (waveguide.k_comp(:, :, 1) .^ 2 ...
%         - (pi / waveguide.a) ^ 2));

    %% SPECTRAL GREEN'S FUNCTIONS
    SGFem(:, :, :, :, media_idx) = dyadic_sgf(lens.er(media_idx), ...
        k(media_idx), k_comp(:, :, :, media_idx), 'E', 'M');

    %% WAVEGUIDE RADIATED FAR-FIELD
    Eff(:, :, :, media_idx) = farfield(k(media_idx), R, sph_grid, ...
        k_comp(:, :, 3, media_idx), SGFem(:, :, :, :, media_idx), M(:, :, :, media_idx));
    Eff_total(:, :, media_idx) = total_field(Eff(:, :, :, media_idx));
end

%% CARTESIAN COORDINATES
sph_coord = ones( [size(sph_grid, 1, 2), 3] ) * 3.21e-3;
sph_coord(:, :, 2 : 3) = sph_grid;
cart_coord = sph2cart_cord(sph_coord);

%% UV COORDINATES
uv_grid = uv_repr(sph_grid);

for media_idx = 1 : 1 : length(lens.er)
    %% PLOT ELECTRIC CURRENT DENSITY
    figure('Position', [250 250 1050 400]);
    subplot(1, 2, 1);
    surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
        norm_magnitude(M(:, :, 1, media_idx), 'dB'), 'LineStyle', 'none');
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
    zlabel('|M_{eq}| / dB');
    subplot(1, 2, 2);
    surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
        norm_magnitude(M(:, :, 1, media_idx), 'dB'), 'LineStyle', 'none');
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
    zlabel('|M_{eq}| / dB');
    sgtitle(['Waveguide M_{eq} @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
        'a = ' num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
        num2str(round(waveguide.b * 1e3, 2)) ' mm'], 'FontSize', ...
        17, 'FontWeight', 'bold');
    saveas(gcf, ['figures\waveguide_Meq_lens_er_' ...
        num2str(lens.er) '.fig']);

    %% PLOT APERTURE ELECTRIC FAR-FIELD
    figure('Position', [250 250 1050 400]);
    subplot(1, 2, 1);
    surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
        norm_magnitude(Eff_total(:, :, media_idx), 'dB'), ...
        'LineStyle', 'none');
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
        norm_magnitude(Eff_total(:, :, media_idx), 'dB'), ...
        'LineStyle', 'none');
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
    sgtitle(['Waveguide E^{FF} @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
        '\epsilon_{r,lens} = ' num2str(lens.er(media_idx)) ', a = ' ...
        num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
        num2str(round(waveguide.b * 1e3, 2)) ' mm'], 'FontSize', ...
        17, 'FontWeight', 'bold');
    saveas(gcf, ['figures\waveguide_pattern_lens_er_' ...
        num2str(lens.er) '.fig']);
end
