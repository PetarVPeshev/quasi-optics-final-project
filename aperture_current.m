% close all;
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
wave.f = 70e9;
waveguide.E10 = 1;
waveguide.er = 1;
waveguide.a = 3.212e-3;
waveguide.b = 1.6e-3;
lens.er = [11.9 4 2];
% lens.max_angle = [40 40 40] * pi / 180;
lens_grid.Nrho = 200;
lens_grid.Nphi = 200;
% n = 4;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;
waveguide.kx1 = pi / waveguide.a;
waveguide.kz = sqrt(wave.k0 ^ 2 - waveguide.kx1 .^ 2);
lens.D = 6 * wave.wavelength;

%% CYLINDRICAL COORDINATE GRID INSIDE LENS
lens_grid.rho = linspace(eps, lens.D / 2, lens_grid.Nrho);
lens_grid.phi = linspace(eps, 2 * pi, lens_grid.Nphi);
lens_grid.cyl_grid = meshgrid_comb(lens_grid.rho, lens_grid.phi);

%% REFRACTIVE INDECIES
waveguide.n = sqrt(waveguide.er);
lens.n = sqrt(lens.er);
lens.Z = zeta ./ lens.n;

%% WAVEGUIDE TE IMPEDANCE
waveguide.ZTE = zeta * wave.k0 ./ waveguide.kz;

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
waveguide_lens.TE_coef = 2 * lens.Z ./ (lens.Z + waveguide.ZTE);
waveguide_lens.TE_T = (waveguide_lens.TE_coef .^ 2) ...
    .* (waveguide.ZTE ./ lens.Z);

%% LENS PARAMETERS
J = zeros( [size(lens_grid.cyl_grid, 1, 2), 3, length(lens.er)] );
J_norm = NaN( [size(J, 1, 2, 3), length(lens.er)] );
for media_idx = 1 : 1 : length(lens.er)
    lens.param(media_idx) = lens_parameters(lens.D, lens.er(media_idx), ...
        lens_grid.cyl_grid);  % , lens.max_angle(media_idx)

    S = sqrt(cos(lens.param(media_idx).theta_t) .* (lens.param(media_idx).eccint ...
        * cos(lens.param(media_idx).sph_grid(:, :, 2)) - 1) ./ ( cos(lens.param(media_idx).theta_i) ...
        .* (lens.param(media_idx).eccint - cos(lens.param(media_idx).sph_grid(:, :, 2))) ));

    [k_comp, k] = wave_vector(lens.er(media_idx), wave.k0, ...
        lens.param(media_idx).sph_grid(:, :, 2 : 3));

    M = zeros( [size(lens.param(media_idx).sph_grid, 1, 2), 3] );
    M(:, :, 1) = 4 * pi * waveguide_lens.TE_T(media_idx) * waveguide.E10 ...
        * waveguide.b * cos(waveguide.a * k_comp(:, :, 1) / 2) ...
        .* sinc(waveguide.b * k_comp(:, :, 2) / 2) ...
        ./ (waveguide.a * (k_comp(:, :, 1) .^ 2 ...
        - (pi / waveguide.a) ^ 2));

    SGFem = dyadic_sgf(lens.er(media_idx), k, k_comp, 'E', 'M');

%     const = ( cos(lens.param(media_idx).sph_grid(:, :, 2)) .^ n ) ./ lens.param(media_idx).sph_grid(:, :, 1);
%     Efeed = zeros( [size(lens.param(media_idx).sph_grid, 1, 2), 3] );
%     Efeed(:, :, 2) = const .* cos(lens.param(media_idx).sph_grid(:, :, 3));
%     Efeed(:, :, 3) = - const .* sin(lens.param(media_idx).sph_grid(:, :, 3));

    Efeed = farfield(k, lens.param(media_idx).sph_grid(:, :, 1), ...
        lens.param(media_idx).sph_grid(:, :, 2 : 3), k_comp, SGFem, M);
    Efeed_total = total_field(Efeed);
%     Efeed = Efeed .* exp(1j * k * lens.param(media_idx).sph_grid(:, :, 1));

    [par_coeff, per_coeff] = transm_coeff(lens.param(media_idx).theta_i, ...
        lens.param(media_idx).theta_t, 1, lens.er(media_idx));

    Ja_const = - 2 * S ./ zeta;
    J(:, :, 1, media_idx) = Ja_const .* par_coeff ...
        .* Efeed(:, :, 2);
    J(:, :, 2, media_idx) = Ja_const .* per_coeff ...
        .* Efeed(:, :, 3);
    J(:, :, :, media_idx) = cyl2cart_vector(J(:, :, :, media_idx), lens_grid.cyl_grid(:, :, 2));

    cyl_grid = zeros( [size(lens_grid.cyl_grid, 1, 2), 3] );
    cyl_grid(:, :, 1 : 2) = lens_grid.cyl_grid;
    cart_grid = cyl2cart_cord(cyl_grid);
    
    J_max_magn = 20 * log10( max(abs(J(:, :, :, media_idx)), [], 'all') );
    J_norm(:, :, 1, media_idx) = 20 * log10( abs(J(:, :, 1, media_idx)) ) - J_max_magn;
    J_norm(:, :, 2, media_idx) = 20 * log10( abs(J(:, :, 2, media_idx)) ) - J_max_magn;
    J_norm(:, :, 3, media_idx) = 20 * log10( abs(J(:, :, 3, media_idx)) ) - J_max_magn;

%     uv_grid = uv_repr(lens.param(media_idx).sph_grid(:, :, 2 : 3));
%     figure('Position', [250 250 1050 400]);
%     subplot(1, 2, 1);
%     surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
%         norm_magnitude(Efeed_total, 'dB'), ...
%         'LineStyle', 'none');
%     grid on;
%     colormap('jet');
%     colorbar;
%     xlim([-1 1]);
%     ylim([-1 1]);
%     caxis([-40 0]);
%     view(0, 90);
%     xticks(-1 : 1 : 1);
%     yticks(-1 : 0.5 : 1);
%     xlabel('U');
%     ylabel('V');
%     zlabel('|E| / dB');
%     subplot(1, 2, 2);
%     surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
%         norm_magnitude(Efeed_total, 'dB'), ...
%         'LineStyle', 'none');
%     grid on;
%     colormap('jet');
%     colorbar;
%     xlim([-1 1]);
%     ylim([-1 1]);
%     caxis([-40 0]);
%     view(-37.5, 30);
%     xticks(-1 : 1 : 1);
%     yticks(-1 : 0.5 : 1);
%     xlabel('U');
%     ylabel('V');
%     zlabel('|E| / dB');
%     sgtitle(['Waveguide E^{FF} @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
%         '\epsilon_{r,lens} = ' num2str(lens.er(media_idx)) ', a = ' ...
%         num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
%         num2str(round(waveguide.b * 1e3, 2)) ' mm'], 'FontSize', ...
%         17, 'FontWeight', 'bold');
end

%% PLOT EQUIVALENT APERTURE ELECTRIC CURRENT DENSITY
for media_idx = 1 : 1 : length(lens.er)
    figure('Position', [250 250 1050 400]);
    subplot(1, 2, 1);
    surface(cart_grid(:, :, 1) * 1e3, cart_grid(:, :, 2) * 1e3, ...
        J_norm(:, :, 1, media_idx), 'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
%     xticks(-5 : 2.5 : 5);
%     yticks(-5 : 2.5 : 5);
%     caxis([-10 0]);
%     caxis([-80 -60]);
    caxis([-80 0]);
    view(0, 90);
    xlabel('x / mm');
    ylabel('y / mm');
    zlabel('|J_{x}| / dB');
    title('|J_{x}| / dB');
    subplot(1, 2, 2);
    surface(cart_grid(:, :, 1) * 1e3, cart_grid(:, :, 2) * 1e3, ...
        J_norm(:, :, 2, media_idx), 'LineStyle', 'none');
    grid on;
    colormap('jet');
    colorbar;
%     xticks(-5 : 2.5 : 5);
%     yticks(-5 : 2.5 : 5);
%     caxis([-45 -25]);
%     caxis([-120 -80]);
    caxis([-50 0]);
    view(0, 90);
    xlabel('x / mm');
    ylabel('y / mm');
    zlabel('|J_{y}| / dB');
    title('|J_{y}| / dB');
    sgtitle(['Untruncated Aperture J_{eq} @ f = ' num2str(wave.f * 1e-9) ...
        ' GHz, D = ' num2str(round(lens.D * 1e2, 2)) ' cm, ' ...
        '\epsilon_{r} = ' num2str(lens.er(media_idx)) ', a = ' ...
        num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
        num2str(round(waveguide.b * 1e3, 2)) ' mm'], 'FontSize', 17, ...
        'FontWeight', 'bold');
    saveas(gcf, ['figures\aperture_Jeq_lens_er_' ...
        num2str(lens.er(media_idx)) '.fig']);
end
