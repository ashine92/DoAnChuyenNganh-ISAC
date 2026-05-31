function [theta_az_hat, theta_el_hat] = estimate_angles(AR_hat, params)
    % estimate_angles: Path-by-path 2D-ESPRIT for near-field tensor array
    
    NRy    = params.NRy;
    NRz    = params.NRz;
    L      = params.L;
    d      = params.d;
    lambda = params.lambda;
    NR     = NRy * NRz;

    Wy = (NRy-1)/2;
    Wz = (NRz-1)/2;
    Vy = Wy + 1;
    Vz = Wz + 1;
    N_ds = Vy * Vz;

    nz_idx = @(nz) nz + (NRz+1)/2;
    ny_idx = @(ny) ny + (NRy+1)/2;
    ant_idx = @(ny,nz) (nz_idx(nz)-1)*NRy + ny_idx(ny);

    vy_range = (0:Wy);
    vz_range = (0:Wz);

    theta_az_hat = zeros(L, 1);
    theta_el_hat = zeros(L, 1);

    for c = 1:L
        a_hat = AR_hat(:, c);
        
        U_s_mat = reshape(a_hat, NRz, NRy);
        U1 = U_s_mat(:, 1:NRy-1); U1 = U1(:);
        U2 = U_s_mat(:, 2:NRy);   U2 = U2(:);
        U3 = U_s_mat(1:NRz-1, :); U3 = U3(:);
        U4 = U_s_mat(2:NRz, :);   U4 = U4(:);
        
        lambda_y = pinv(U1) * U2;
        lambda_z = pinv(U3) * U4;
        
        % The phase of nearfield array response is +2*pi/lambda * d * (...)
        % So phase diff is +pi * sin_el * sin_az and +pi * cos_el
        el = acos(max(-1, min(1, angle(lambda_z) / pi)));
        sin_el = max(eps, sin(el));
        az = asin(max(-1, min(1, angle(lambda_y) / pi ./ sin_el)));
        
        theta_el_hat(c) = real(el);
        theta_az_hat(c) = real(az);
    end
end
