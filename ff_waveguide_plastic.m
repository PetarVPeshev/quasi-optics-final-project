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
lens.er = 2;
% Grid parameters
Ntheta = 800;
Nphi = 3200;
% Far-field parameters
R = 1;
% Plane plots
planes = [0 90];

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
lens.n = sqrt(lens.er);
lens.Z = zeta / lens.n;

%% TRANSMISSION COEFFICIENT @ WAVEGUIDE-LENS INTERFACE
waveguide_lens.TE_coef = 2 * lens.Z ./ (lens.Z + waveguide.ZTE);

%% COORDINATE GRID
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = linspace(0, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

%% WAVE VECTOR COMPONENTS
[k_comp, k] = wave_vector(lens.er, wave.k0, sph_grid);

%% EQUIVALENT MAGNETIC CURRENT AND RADIATED FAR-FIELD
% Waveguide radiated field and equivalent magnetic current
[E, M] = waveguide_feed(waveguide, waveguide_lens.TE_coef, k, k_comp, ...
    R, sph_grid);
% Waveguide total radiated magnetic field
E_total = total_field(E);

%% DIRECTIVITY
dir = directivity(lens.er, E, sph_grid, R);
dir_broadside = dir(1, 1);

%% UV COORDINATES
uv_grid = uv_repr(sph_grid);

%% PLOT MAGNETIC CURRENT DENSITY
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
    norm_magnitude(M(:, :, 1), 'dB'), 'LineStyle', 'none');
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
    norm_magnitude(M(:, :, 1), 'dB'), 'LineStyle', 'none');
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
    '\epsilon_{r} = ' num2str(lens.er) ', a = ' ...
    num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
    num2str(round(waveguide.b * 1e3, 2)) ' mm'], 'FontSize', ...
    17, 'FontWeight', 'bold');
saveas(gcf, 'figures\wg_Meq_plastic.fig');

%% PLOT APERTURE ELECTRIC FAR-FIELD
figure('Position', [250 250 1050 400]);
subplot(1, 2, 1);
surface(uv_grid(:, :, 1), uv_grid(:, :, 2), ...
    norm_magnitude(E_total, 'dB'), 'LineStyle', 'none');
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
    norm_magnitude(E_total, 'dB'), 'LineStyle', 'none');
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
    '\epsilon_{r} = ' num2str(lens.er) ', a = ' ...
    num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
    num2str(round(waveguide.b * 1e3, 2)) ' mm'], 'FontSize', ...
    17, 'FontWeight', 'bold');
saveas(gcf, 'figures\wg_pattern_plastic.fig');

%% PLOT E-PLANE
E_norm = norm_magnitude(E_total, 'dB');
theta_plot = NaN(1, 2 * length(theta));
theta_plot(1 : length(theta)) = - fliplr(theta) * 180 / pi;
theta_plot(length(theta) + 1 : end) = theta * 180 / pi;
line_styles = ["-", "--"];
figure('Position', [250 250 725 400]);
for idx = 1 : 1 : length(planes)
    plane_idx_1 = find(round(phi * 180 / pi, 0) == planes(idx), 1);
    plane_idx_2 = find(round(phi * 180 / pi, 0) == planes(idx) + 180, 1);
    plane_field = NaN(1, length(theta) * 2);
    plane_field(1 : length(theta)) = fliplr(E_norm(plane_idx_2, :));
    plane_field(length(theta) + 1 : end) = E_norm(plane_idx_1, :);
    plot(theta_plot, plane_field, [line_styles(idx)], 'LineWidth', 2.0, ...
        'DisplayName', ['E_{WG}^{FF}, \phi = ' num2str(planes(idx)) ...
        ' deg and ' num2str(planes(idx) + 180) ' deg']);
    hold on;
end
hold off;
grid on;
ylim([-40 0]);
xlim([-90 90]);
xticks(-90 : 15 : 90);
legend show;
legend('location', 'bestoutside');
xlabel('\theta / deg');
ylabel('|E_{WG}^{FF}| / dB');
title(['Waveguide E^{FF} @ f = ' num2str(wave.f * 1e-9) ' GHz, ' ...
    '\epsilon_{r} = ' num2str(lens.er) ', a = ' ...
    num2str(round(waveguide.a * 1e3, 2)) ' mm, and b = ' ...
    num2str(round(waveguide.b * 1e3, 2)) ' mm']);
saveas(gcf, 'figures\wg_pattern_planes_plastic.fig');

%% SAVE WORKSPACE
wg_field_plastic.sph_grid = sph_grid;
wg_field_plastic.theta = theta;
wg_field_plastic.phi = phi;
wg_field_plastic.waveguide = waveguide;
wg_field_plastic.E = E;
wg_field_plastic.M = M;
wg_field_plastic.E_total = E_total;
save('results\wg_field_plastic.mat', 'wg_field_plastic');
