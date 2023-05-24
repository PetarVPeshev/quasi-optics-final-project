function [J, lens, varargout] = lens_current(lens, E_inc, theta_inc, ...
    theta_tr, medium_er, cyl_grid, sph_grid, varargin)
%LENS_CURRENT Summary of this function goes here
%   Detailed explanation goes here
    zeta = 376.730313668 / sqrt(medium_er);
    enable_ft = false;
    if ~isempty(varargin)
        if strcmp(varargin{1}, 'FT')
            enable_ft = true;
            k_comp = varargin{2};
        else
            error('Error. Invalid arguments.');
        end
    elseif length(varargin) > 3
        error('Error. Invalid number of arguments.');
    end
    
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

    if enable_ft == true
        rho = cyl_grid(:, :, 1);
        drho = rho(1, 2) - rho(1, 1);
        phi = cyl_grid(:, :, 2);
        dphi = phi(2, 1) - phi(1, 1);

        num_pages = floor( memory().MaxPossibleArrayBytes * 0.8 / (96 * size(J, 1) * size(J, 2)) );
        while rem(size(k_comp, 1) * size(k_comp, 2), num_pages) ~= 0
            num_pages = num_pages - 1;
        end
        num_k_elements = size(k_comp, 1) * size(k_comp, 2);
        num_iterations = num_k_elements / num_pages;
        
        Jft = zeros( [size(k_comp, 1, 2), 3] );
        for coord_idx = 1 : 1 : 2
            Jcoord_idx = NaN(num_k_elements, 1);
            for iteration_idx = 1 : 1 : num_iterations
                start_idx = (iteration_idx - 1) * num_pages + 1;
                end_idx = iteration_idx * num_pages;

                kx = k_comp(start_idx : end_idx);
                ky = k_comp(num_k_elements + start_idx : num_k_elements + end_idx);
                kx = permute(repmat(kx, [1 1 size(J, 1) size(J, 2)]), [3 4 2 1]);
                ky = permute(repmat(ky, [1 1 size(J, 1) size(J, 2)]), [3 4 2 1]);
                
                Jiteration = sum( sum( J(:, :, coord_idx) .* exp( 1j * kx .* rho .* cos(phi) ) .* exp( 1j * ky .* rho .* sin(phi) ) .* rho ) ) * drho * dphi;
                Jcoord_idx(start_idx : end_idx) = permute(Jiteration(1, 1, :), [3 1 2]);
            end
            Jft(:, :, coord_idx) = reshape(Jcoord_idx, [size(k_comp, 1) size(k_comp, 2)]);
        end

        varargout{1} = Jft;
    end
end

