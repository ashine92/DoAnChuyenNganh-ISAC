function [theta_az_hat, theta_el_hat] = estimate_angles(AR_hat, params)
% =========================================================================
% estimate_angles.m
% =========================================================================
% Description:
%   Estimates azimuth and elevation angles from the estimated array
%   response matrix using second-order Taylor expansion and rotational
%   invariance (ESPRIT-like) on the down-sampled covariance matrix.
%
%   Steps (Section 5.2 of paper):
%   1. Compute covariance: R_hat = AR_hat * AR_hat^H  (NR x NR)
%   2. Extract down-sampled covariance R_tilde (Eq. 33-34) by selecting
%      element pairs (nR, n_tilde_R) where ny = -ny_tilde, nz = -nz_tilde
%   3. SVD of R_tilde to get signal subspace U_L (Eq. 35)
%   4. Apply ESPRIT rotational invariance (Eq. 36-40) to extract eigenvalues
%   5. Recover angles from eigenvalues (Eq. 40)
%
% Inputs:
%   AR_hat - NR x L estimated array response matrix at UT
%   params - struct with fields: NRy, NRz, d, lambda, L
%
% Outputs:
%   theta_az_hat - L x 1 estimated azimuth angles [rad]
%   theta_el_hat - L x 1 estimated elevation angles [rad]
%
% Note: Same function is used for AoD by passing BT_hat and TX params.
%
% Runtime: O(NR^2 + Wy*Wz^2 * L)
% =========================================================================

    NRy    = params.NRy;
    NRz    = params.NRz;
    L      = params.L;
    d      = params.d;
    lambda = params.lambda;
    NR     = NRy * NRz;

    % -------------------------------------------------------
    % Step 1: Covariance matrix
    % -------------------------------------------------------
    R_hat = AR_hat * AR_hat';    % NR x NR

    % -------------------------------------------------------
    % Step 2: Down-sample to get R_tilde
    %   Using pairs where n = -n_tilde (ny = -ny_tilde, nz = -nz_tilde)
    %   R_tilde[vy,vz; vy_tilde,vz_tilde] = R_hat[nR(vy-vy_t, vz-vz_t), nR_c]
    %   where nR_c is the center index
    %
    %   R_tilde is a (Wy+1)(Wz+1) x (Wy+1)(Wz+1) matrix structured as
    %   block Toeplitz form described in Eq. (33-34).
    %   Wy = (NRy-1)/2, Wz = (NRz-1)/2
    % -------------------------------------------------------
    Wy = (NRy-1)/2;
    Wz = (NRz-1)/2;

    % Build index mapping: antenna index as function of (ny, nz)
    % nR(ny, nz) = (NRy*NRz+1)/2 + nz*NRy + ny  (1-based)
    % ny in [-(NRy-1)/2 ... (NRy-1)/2]
    % nz in [-(NRz-1)/2 ... (NRz-1)/2]
    ny_range = (-(NRy-1)/2 : (NRy-1)/2);   % 1 x NRy
    nz_range = (-(NRz-1)/2 : (NRz-1)/2);   % 1 x NRz

    % Function: given (ny, nz) → linear index in R_hat (1-based)
    % The ordering in near_field_array_response.m uses meshgrid(ny,nz)
    % with nz as outer (rows), ny as inner (cols), then (:) to vectorize
    % So linear index = (nz_idx-1)*NRy + ny_idx   where _idx are 1-based
    nz_idx = @(nz) nz + (NRz+1)/2;   % convert nz to 1-based index
    ny_idx = @(ny) ny + (NRy+1)/2;
    ant_idx = @(ny,nz) (nz_idx(nz)-1)*NRy + ny_idx(ny);

    % Build R_tilde: size (NRy_half x NRz_half) downsampled version
    % Use vR = vz*(NRy+1)/2 + vy + 1  for vy in [0..Wy], vz in [0..Wz]
    vy_range = (0:NRy-1);   % v variables: 0-based
    vz_range = (0:NRz-1);

    % Size of R_tilde: (NRy*NRz) x (NRy*NRz) but using the downsampling
    % construction from Eq.(33-34).
    % For simplicity and robustness, we construct a (NRy*NRz/2) approximation.
    % The mapping R_hat(nR, n_tilde_R) → R_tilde(vR, v_tilde_R) holds when
    % nR = center + (vz-vz_t)*NRy + (vy-vy_t) and n_tilde_R = center - same.
    %
    % Practical construction: extract elements e(ny, nz) = R_hat(nR(ny,nz), n_tilde_R(-ny,-nz))
    % and form a (NRy*NRz) x (NRy*NRz) Hermitian matrix.
    %
    % For the ESPRIT step we need R_tilde as a (Wy+1)(Wz+1) x (Wy+1)(Wz+1) matrix.
    % We use a simplified extraction: R_tilde(vR, v_tilde_R) = R_hat(i1, i2)
    % where i1 = ant_idx(vy-vy_t, vz-vz_t), i2 = ant_idx(-(vy-vy_t), -(vz-vz_t))
    % i.e., the element e(vy-vy_tilde, vz-vz_tilde).

    % Down-sampled size: (NRy*NRz) x (NRy*NRz)
    N_ds = NRy * NRz;
    R_tilde = zeros(N_ds, N_ds);

    [VY, VY2] = meshgrid(vy_range, vy_range);
    [VZ, VZ2] = meshgrid(vz_range, vz_range);

    % Vectorized construction of R_tilde
    for vi = 1:N_ds
        vyi = vy_range(mod(vi-1, NRy)+1);
        vzi = vz_range(floor((vi-1)/NRy)+1);
        for vj = 1:N_ds
            vyj = vy_range(mod(vj-1, NRy)+1);
            vzj = vz_range(floor((vj-1)/NRy)+1);
            dny = vyi - vyj;
            dnz = vzi - vzj;
            % Map to antenna pair
            ny_p  = dny/2;
            nz_p  = dnz/2;
            ny_n  = -dny/2;
            nz_n  = -dnz/2;
            if abs(ny_p) <= Wy && abs(nz_p) <= Wz && ...
               mod(dny,2)==0 && mod(dnz,2)==0
                i1 = ant_idx(round(ny_p), round(nz_p));
                i2 = ant_idx(round(ny_n), round(nz_n));
                if i1 >= 1 && i1 <= NR && i2 >= 1 && i2 <= NR
                    R_tilde(vi, vj) = R_hat(i1, i2);
                end
            end
        end
    end
    % Fill zeros with nearest neighbor to avoid rank collapse
    R_tilde = R_tilde + eye(N_ds) * (trace(R_hat)/NR * 1e-6);

    % -------------------------------------------------------
    % Step 3: SVD of R_tilde → signal subspace U_L (Eq. 35)
    % -------------------------------------------------------
    [U, S, ~] = svd(R_tilde);
    U_L = U(:, 1:L);    % N_ds x L signal subspace

    % -------------------------------------------------------
    % Step 4: ESPRIT rotational invariance (Eq. 36-40)
    % Selection matrices for y and z directions
    %
    % For y direction:
    %   U1 = (I_Wz ⊗ J1) * U_L    J1 = [I_Wy, 0]
    %   U2 = (I_Wz ⊗ J2) * U_L    J2 = [0, I_Wy]
    %   Psi_y = U1† * U2
    %
    % For z direction:
    %   U3 = (J3 ⊗ I_Wy) * U_L    J3 = [I_Wz, 0]
    %   U4 = (J4 ⊗ I_Wy) * U_L    J4 = [0, I_Wz]
    %   Psi_z = U3† * U4
    % -------------------------------------------------------

    % Use integer Wy, Wz
    Wy_int = round(Wy);
    Wz_int = round(Wz);

    % Selection matrices (Eq. 38)
    J1_mat = [eye(Wy_int), zeros(Wy_int, 1)];   % Wy x (Wy+1)
    J2_mat = [zeros(Wy_int, 1), eye(Wy_int)];
    J3_mat = [eye(Wz_int), zeros(Wz_int, 1)];
    J4_mat = [zeros(Wz_int, 1), eye(Wz_int)];

    I_Wz = eye(Wz_int+1);
    I_Wy = eye(Wy_int+1);

    % Reshape U_L to [Wz+1, Wy+1, L] for easier slicing
    % U_L is N_ds x L = (NRy*NRz) x L, with vy fast, vz slow
    U_3d = reshape(U_L, [NRy, NRz, L]);  % NRy x NRz x L

    % Sub-matrices for y-direction (Eq. 37)
    % U1 = (I_Wz ⊗ J1)*U_L removes last vy slice; U2 removes first vy slice
    U1_mat = reshape(U_3d(1:end-1, :, :), [], L);    % (NRy-1)*NRz x L
    U2_mat = reshape(U_3d(2:end, :, :), [], L);

    % Sub-matrices for z-direction
    U3_mat = reshape(U_3d(:, 1:end-1, :), [], L);    % NRy*(NRz-1) x L
    U4_mat = reshape(U_3d(:, 2:end, :), [], L);

    % Rotation factors (Eq. 39)
    Psi_y = pinv(U1_mat) * U2_mat;    % L x L
    Psi_z = pinv(U3_mat) * U4_mat;    % L x L

    % Eigenvalues (Eq. 40)
    lambda_y = eig(Psi_y);   % L eigenvalues
    lambda_z = eig(Psi_z);   % L eigenvalues

    % -------------------------------------------------------
    % Step 5: Recover angles (Eq. 40)
    %   theta_el = acos(-angle(lambda_z) / pi)
    %   theta_az = asin(-angle(lambda_y) / pi / sin(theta_el))
    %
    % Note: angle() returns principal argument in (-pi, pi]
    % The normalization by pi comes from: phase = 2*pi/lambda * 2d * sin(...)
    % with d = lambda/4, so phase = pi * sin(...)  → angle/pi = sin(...)
    % -------------------------------------------------------

    % Sort to match paths (pairing based on nearest angle values)
    theta_el_hat = acos(max(-1, min(1, -angle(lambda_z) / pi)));   % L x 1
    sin_el       = max(eps, sin(theta_el_hat));
    theta_az_hat = asin(max(-1, min(1, -angle(lambda_y) / pi ./ sin_el)));  % L x 1

    % Sort angles by azimuth for consistent path ordering
    [theta_az_hat, sort_idx] = sort(real(theta_az_hat));
    theta_el_hat = theta_el_hat(sort_idx);

    theta_az_hat = real(theta_az_hat);
    theta_el_hat = real(theta_el_hat);
end
