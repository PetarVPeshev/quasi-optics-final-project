function [E, M] = waveguide_feed(waveguide, te_coeff, k, k_comp, r, ...
    sph_grid, varargin)
%WAVEGUIDE_FEED Summary of this function goes here
%   Detailed explanation goes here

    a = waveguide.a;
    b = waveguide.b;
    E10 = waveguide.E10;
    
    kx = k_comp(:, :, 1);
    ky = k_comp(:, :, 2);
    kz = k_comp(:, :, 3);
    
    theta = sph_grid(:, :, 1);
    phi = sph_grid(:, :, 2);
    
    Mx = 4 * pi * te_coeff * E10 * (b / a) * cos(a * kx / 2) ...
        .* sinc(b * ky / (2 * pi)) ./ (kx .^ 2 - (pi / a) ^ 2);
    
    % xx-component of SGFem is equal to zero
    G_const = - 1 ./ (2 * kz);
    Gyx = - 1j * G_const .* kz;
    Gzx = 1j * G_const .* ky;
    
    E_const = 1j * kz .* exp(-1j * k * r) ./ (2 * pi * r);
    if ~isempty(varargin)
        if strcmp(varargin{1}, 'NeglectPhase')
            E_const = 1j * kz ./ (2 * pi * r);
        end
    end
    Ey = E_const .* Gyx .* Mx;
    Ez = E_const .* Gzx .* Mx;
    
    E = zeros( [size(sph_grid, 1, 2), 3] );
    % Theta component
    E(:, :, 2) = Ey .* cos(theta) .* sin(phi) - Ez .* sin(theta);
    % Phi component
    E(:, :, 3) = Ey .* cos(phi);
    
    M = zeros( [size(sph_grid, 1, 2), 3] );
    M(:, :, 1) = Mx;
end

