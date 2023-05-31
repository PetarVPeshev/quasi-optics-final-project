function [E, M] = waveguide_feed(waveguide, te_coeff, k, k_comp, r, ...
    sph_grid, varargin)
%WAVEGUIDE_FEED Summary of this function goes here
%   Detailed explanation goes here

    a = waveguide.a;
    b = waveguide.b;
    E10 = waveguide.E10;
    
    kx = k_comp(:, :, 1);
    ky = k_comp(:, :, 2);
    
    theta = sph_grid(:, :, 1);
    phi = sph_grid(:, :, 2);
    
    Mx = 4 * pi * te_coeff * E10 * (b / a) * cos(a * kx / 2) ...
        .* sinc(b * ky / (2 * pi)) ./ (kx .^ 2 - (pi / a) ^ 2);
    
    E_const = 1j * k * Mx .* exp(-1j * k * r) ./ (4 * pi * r);
    if ~isempty(varargin)
        if strcmp(varargin{1}, 'NeglectPhase')
            E_const = 1j * k * Mx ./ (4 * pi * r);
        end
    end
    E = zeros( [size(sph_grid, 1, 2), 3] );
    % Theta component
    E(:, :, 2) = E_const .* sin(phi);
    % Phi component
    E(:, :, 3) = E_const .* cos(theta) .* cos(phi);
    
    M = zeros( [size(sph_grid, 1, 2), 3] );
    M(:, :, 1) = Mx;
end

