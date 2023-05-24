function [J, lens] = lens_current(lens, E_inc, theta_inc, theta_tr, medium_er, cyl_grid, sph_grid)
%LENS_CURRENT Summary of this function goes here
%   Detailed explanation goes here
    zeta = 376.730313668 / sqrt(medium_er);
    
    theta = sph_grid(:, :, 1);
    
    [TM_coef, TE_coef] = transm_coeff(theta_inc, theta_tr, medium_er, lens.er);
    
    S = sqrt(cos(theta_tr) .* (lens.e * cos(theta) - 1) ...
        ./ ( cos(theta_inc) .* (lens.e - cos(theta)) ));
    
    J_const = - 2 * S / zeta;
    J = zeros( [size(J_const, 1, 2), 3] );
    J(:, :, 1) = J_const .* TM_coef .* E_inc(:, :, 2);
    J(:, :, 2) = J_const .* TE_coef .* E_inc(:, :, 3);
    J = cyl2cart_vector(J, cyl_grid);
    
    lens.TM_coef = TM_coef;
    lens.TE_coef = TE_coef;
end

