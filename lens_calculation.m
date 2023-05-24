function [lens, cyl_grid, sph_grid, theta_inc, theta_tr] = ...
    lens_calculation(lens, cyl_grid, varargin)
%LENS_CALCULATION Summary of this function goes here
%   Detailed explanation goes here
    rho = cyl_grid(:, :, 1);
    phi = cyl_grid(:, :, 2);

    lens.R = lens.D / 2;                    % Radius
    lens.e = 1 / sqrt(lens.er);             % Eccintercity
    lens.theta_crit = asin(lens.e);         % Reflection critical angle
    
    % Maximum angle, equals critical angle if not specified
    if isempty(varargin)
        lens.theta_max = pi / 2 - lens.theta_crit - 0.01 * pi / 180;
    elseif length(varargin) == 1
        lens.theta_max = varargin{1};
    else
        error('Error. Invalid number of arguments.');
    end
    
    % Minimum R, always @ maximum angle
    lens.r_min = lens.R / sin(lens.theta_max);
    
    % Semi-major axis
    lens.a = lens.r_min * (1 - lens.e * cos(lens.theta_max)) / (1 - lens.e ^ 2);

    % Foci distance
    lens.c = lens.a * lens.e;

    % Semi-minor axis
    lens.b = sqrt(lens.a ^ 2 - lens.c ^ 2);
    
    % Cylindrical coordinates z-axis
    z = lens.a * sqrt(1 - (rho / lens.b) .^ 2) + lens.c;

    cyl_grid(:, :, 3) = z;

    % Spherical coordinates transformation
    theta = atan(rho ./ z);
    r = lens.a * (1 - lens.e ^ 2) ./ (1 - lens.e * cos(theta));

    sph_grid = NaN( [size(cyl_grid, 1, 2), 3] );
    sph_grid(:, :, 1) = r;
    sph_grid(:, :, 2) = theta;
    sph_grid(:, :, 3) = phi;

    % Incident & transmission angle @ lens dielectric-free space interface
    theta_inc = acos( (1 - lens.e * cos(theta)) ./ sqrt(1 + lens.e ^ 2 - 2 * lens.e * cos(theta)) );
    theta_tr = acos( (cos(theta) - lens.e) ./ sqrt(1 + lens.e ^ 2 - 2 * lens.e * cos(theta)) );
%     theta_tr = asin( sqrt(lens.er) * sin(theta_inc) );
end

