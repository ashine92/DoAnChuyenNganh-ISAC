function [CRB_params, CRB_pR, CRB_pl] = compute_crb(params, AR, BT, C, pl_all, SNR_dB)
% =========================================================================
% compute_crb.m  —  Cramér-Rao Bound for channel parameters and positions
% =========================================================================
% Full vectorized FIM:
%   Ω(p,q)(l1,l2) = (2/σ²) Re{ vec(∂S/∂ϖ_p,l1)^H * vec(∂S/∂ϖ_q,l2) }
%
% where S = Σ_l a_l ⊗ b_l ⊗ c_l is the noiseless tensor signal.
%
% σ² is the ACTUAL noise variance per element (matching construct_tensor).
% =========================================================================

    L      = params.L;
    NR     = params.NR;   NT = params.NT;
    NRy    = params.NRy;  NRz = params.NRz;
    NTy    = params.NTy;  NTz = params.NTz;
    K      = params.K;    K_bar = params.K_bar;
    F_size = params.F;    T_size = params.T;
    fs     = params.fs;   d = params.d;
    lambda = params.lambda; c_speed = params.c;
    pT     = params.pT;   pR = params.pR;

    % --- Noise variance: MUST match construct_tensor.m ---
    % In construct_tensor: sig_pwr = ||Y_clean||^2/numel(Y_clean)
    %                      noise_var = sig_pwr / SNR_lin
    % For CRB, we compute sig_pwr from the true factors.
    SNR_lin = 10^(SNR_dB/10);

    % Subcarrier indices
    k_idx = round(linspace(1, K_bar, K))';

    % DFT beamforming matrices (same as construct_tensor)
    if F_size <= NR
        W = dftmtx(NR)/sqrt(NR);
        W = W(:, 1:F_size);
    else
        W = dftmtx(F_size)/sqrt(F_size);
        W = W(1:NR, :);
    end

    if T_size <= NT
        F_mat = dftmtx(NT)/sqrt(NT);
        F_mat = F_mat(:, 1:T_size);
    else
        F_mat = dftmtx(T_size)/sqrt(T_size);
        F_mat = F_mat(1:NT, :);
    end

    % Effective factor matrices in beamspace
    A = W' * AR;     % F×L
    B = F_mat' * BT; % T×L

    % Compute actual signal power (matching construct_tensor.m)
    % Y_clean = Σ_l a_l ⊗ b_l ⊗ c_l, so ||Y_clean||²_F = ...
    sig_pwr_total = 0;
    for k = 1:K
        Yk = zeros(F_size, T_size);
        for l = 1:L
            Yk = Yk + A(:,l) * B(:,l)' * C(k,l);
        end
        sig_pwr_total = sig_pwr_total + norm(Yk, 'fro')^2;
    end
    sig_pwr_per_elem = sig_pwr_total / (F_size * T_size * K);
    sigma2 = sig_pwr_per_elem / SNR_lin;

    % Antenna index vectors
    [NY_R, NZ_R] = meshgrid(-(NRy-1)/2:(NRy-1)/2, -(NRz-1)/2:(NRz-1)/2);
    ny_R = NY_R(:);  nz_R = NZ_R(:);

    [NY_T, NZ_T] = meshgrid(-(NTy-1)/2:(NTy-1)/2, -(NTz-1)/2:(NTz-1)/2);
    ny_T = NY_T(:);  nz_T = NZ_T(:);

    % Recover true geometry
    az_R = zeros(L,1); el_R = zeros(L,1);
    az_T = zeros(L,1); el_T = zeros(L,1);
    tau  = zeros(L,1);
    dcR  = zeros(L,1); dcT = zeros(L,1);

    for l = 1:L
        pl      = pl_all(:,l);
        dcR(l)  = norm(pl - pR);
        dcT(l)  = norm(pl - pT);
        el_R(l) = acos(max(-1,min(1,(pl(3)-pR(3))/dcR(l))));
        el_T(l) = acos(max(-1,min(1,(pl(3)-pT(3))/dcT(l))));

        sin_el_R_l = sin(el_R(l));
        if abs(sin_el_R_l) > 1e-10
            sin_az_R_l = (pl(2)-pR(2)) / (dcR(l) * sin_el_R_l);
            sin_az_R_l = max(-1, min(1, sin_az_R_l));
            az_R(l) = atan2(sin_az_R_l, sqrt(max(0, 1-sin_az_R_l^2)));
        end
        sin_el_T_l = sin(el_T(l));
        if abs(sin_el_T_l) > 1e-10
            sin_az_T_l = (pl(2)-pT(2)) / (dcT(l) * sin_el_T_l);
            sin_az_T_l = max(-1, min(1, sin_az_T_l));
            az_T(l) = atan2(sin_az_T_l, sqrt(max(0, 1-sin_az_T_l^2)));
        end
        tau(l)  = (dcR(l)+dcT(l))/c_speed;
    end

    % =====================================================================
    % Build FULL FIM using vectorized signal derivatives
    % =====================================================================
    % Signal model: Y(:,:,k) = A * diag(C(k,:)) * B'
    % vec(Y) = Σ_l (conj(b_l) ⊗ a_l) * c_l(k) stacked over k
    %
    % Derivative wrt parameter ϖ_{p,l}:
    %   ∂vec(Y)/∂ϖ_{p,l} = stacked over k of:
    %     c_l(k) * (conj(b_l) ⊗ ∂a_l/∂ϖ)  for AoA params
    %     c_l(k) * (conj(∂b_l/∂ϖ) ⊗ a_l)  for AoD params  
    %     ∂c_l(k)/∂τ * (conj(b_l) ⊗ a_l)   for ToA
    %
    % FIM:  Ω(i,j) = (2/σ²) Re{ d_i^H * d_j }
    % where d_i = vec(∂Y/∂ϖ_i) over all F,T,K
    % =====================================================================

    % Precompute derivatives of steering vectors
    dAR_az = zeros(NR, L);  dAR_el = zeros(NR, L);
    dBT_az = zeros(NT, L);  dBT_el = zeros(NT, L);

    for l = 1:L
        pl = pl_all(:,l);

        % Distances from SP to each antenna
        ant_R = pR + [zeros(1,NR); ny_R'*d; nz_R'*d];
        dnR   = sqrt(sum((ant_R - pl).^2, 1))';

        ant_T = pT + [zeros(1,NT); ny_T'*d; nz_T'*d];
        dnT   = sqrt(sum((ant_T - pl).^2, 1))';

        % ∂a_R/∂az_R: derivative of near-field steering vector
        % a_R(n) = exp(-j*2π/λ*(d_n - d_c))
        % ∂d_n/∂az = d * ny * cos(az)*sin(el) * d_c / d_n  (approx)
        pf_az = (-1j*2*pi/lambda) * d * cos(az_R(l))*sin(el_R(l));
        dAR_az(:,l) = pf_az * ny_R .* AR(:,l) .* (dcR(l)./(dnR+eps));

        pf_el = (-1j*2*pi/lambda) * d;
        dAR_el(:,l) = pf_el * (sin(az_R(l))*cos(el_R(l))*ny_R ...
                      - sin(el_R(l))*nz_R) .* AR(:,l) .* (dcR(l)./(dnR+eps));

        pf_az_T = (-1j*2*pi/lambda) * d * cos(az_T(l))*sin(el_T(l));
        dBT_az(:,l) = pf_az_T * ny_T .* BT(:,l) .* (dcT(l)./(dnT+eps));

        pf_el_T = (-1j*2*pi/lambda) * d;
        dBT_el(:,l) = pf_el_T * (sin(az_T(l))*cos(el_T(l))*ny_T ...
                      - sin(el_T(l))*nz_T) .* BT(:,l) .* (dcT(l)./(dnT+eps));
    end

    % Transform to beamspace
    dA_az  = W' * dAR_az;   % F×L
    dA_el  = W' * dAR_el;
    dB_az  = F_mat' * dBT_az; % T×L
    dB_el  = F_mat' * dBT_el;

    % dC/dtau
    dC_tau = zeros(K, L);
    for l = 1:L
        dC_tau(:,l) = -1j*2*pi*fs/K_bar * k_idx .* C(:,l);
    end

    % =====================================================================
    % Build full FIM using vectorized derivatives
    % For each parameter pair (p_l1, q_l2), compute:
    %   Ω(p_l1, q_l2) = (2/σ²) * Re{ Σ_k d_{p,l1}(k)^H * d_{q,l2}(k) }
    % where d_{az_R,l}(k) = c_l(k) * kron(conj(B(:,l)), dA_az(:,l))
    %       d_{el_R,l}(k) = c_l(k) * kron(conj(B(:,l)), dA_el(:,l))
    %       d_{az_T,l}(k) = c_l(k) * kron(conj(dB_az(:,l)), A(:,l))
    %       d_{el_T,l}(k) = c_l(k) * kron(conj(dB_el(:,l)), A(:,l))
    %       d_{tau,l}(k)  = dC_tau(k,l) * kron(conj(B(:,l)), A(:,l))
    %
    % Since kron(a,b)^H * kron(c,d) = (a^H*c) * (b^H*d), we can factorize:
    %   Σ_k d_i(k)^H * d_j(k) = [spatial_part] * [subcarrier_part]
    % =====================================================================

    Omega = zeros(5*L, 5*L);

    % Precompute subcarrier inner products: Σ_k c*(k,l1)*c(k,l2)
    CC = C' * C;             % L×L: CC(l1,l2) = Σ_k c*(k,l1)*c(k,l2)
    CdC = C' * dC_tau;       % L×L: CdC(l1,l2) = Σ_k c*(k,l1)*dC(k,l2)
    dCdC = dC_tau' * dC_tau; % L×L

    % Spatial inner products
    AA = A' * A;               % L×L
    BB = B' * B;               % L×L (note: conj(B)^H*conj(B) = B'*B)
    
    dAaz_A  = dA_az' * A;     dAaz_dAaz = dA_az' * dA_az;
    dAel_A  = dA_el' * A;     dAel_dAel = dA_el' * dA_el;
    dAaz_dAel = dA_az' * dA_el;
    
    dBaz_B  = dB_az' * B;     dBaz_dBaz = dB_az' * dB_az;
    dBel_B  = dB_el' * B;     dBel_dBel = dB_el' * dB_el;
    dBaz_dBel = dB_az' * dB_el;
    
    dAaz_B = dA_az;  % just the matrix itself for cross terms

    % Fill FIM blocks
    for l1 = 1:L
        for l2 = 1:L
            % Subcarrier factor for same-derivative-mode terms
            cc = CC(l1,l2);  % Σ_k c*(l1,k)*c(l2,k)

            % --- Block (az_R, az_R) ---
            % d = c_l(k) * kron(conj(B(:,l)), dA_az(:,l))
            % d1^H * d2 = Σ_k c*(l1)*c(l2) * (B(:,l1)^T*conj(B(:,l2))) * (dA_az(:,l1)^H*dA_az(:,l2))
            val = cc * conj(BB(l1,l2)) * dAaz_dAaz(l1,l2);
            Omega(0*L+l1, 0*L+l2) = (2/sigma2) * real(val);

            % --- Block (el_R, el_R) ---
            val = cc * conj(BB(l1,l2)) * dAel_dAel(l1,l2);
            Omega(1*L+l1, 1*L+l2) = (2/sigma2) * real(val);

            % --- Block (az_T, az_T) ---
            val = cc * conj(dBaz_dBaz(l1,l2)) * AA(l1,l2);
            Omega(2*L+l1, 2*L+l2) = (2/sigma2) * real(val);

            % --- Block (el_T, el_T) ---
            val = cc * conj(dBel_dBel(l1,l2)) * AA(l1,l2);
            Omega(3*L+l1, 3*L+l2) = (2/sigma2) * real(val);

            % --- Block (tau, tau) ---
            val = dCdC(l1,l2) * conj(BB(l1,l2)) * AA(l1,l2);
            Omega(4*L+l1, 4*L+l2) = (2/sigma2) * real(val);

            % --- Cross: (az_R, el_R) ---
            val = cc * conj(BB(l1,l2)) * dAaz_dAel(l1,l2);
            Omega(0*L+l1, 1*L+l2) = (2/sigma2) * real(val);
            Omega(1*L+l2, 0*L+l1) = (2/sigma2) * real(val);

            % --- Cross: (az_T, el_T) ---
            val = cc * conj(dBaz_dBel(l1,l2)) * AA(l1,l2);
            Omega(2*L+l1, 3*L+l2) = (2/sigma2) * real(val);
            Omega(3*L+l2, 2*L+l1) = (2/sigma2) * real(val);

            % --- Cross: (az_R, tau) ---
            val = CdC(l1,l2) * conj(BB(l1,l2)) * (dA_az(:,l1)' * A(:,l2));
            Omega(0*L+l1, 4*L+l2) = (2/sigma2) * real(val);
            Omega(4*L+l2, 0*L+l1) = (2/sigma2) * real(val);

            % --- Cross: (el_R, tau) ---
            val = CdC(l1,l2) * conj(BB(l1,l2)) * (dA_el(:,l1)' * A(:,l2));
            Omega(1*L+l1, 4*L+l2) = (2/sigma2) * real(val);
            Omega(4*L+l2, 1*L+l1) = (2/sigma2) * real(val);

            % --- Cross: (az_T, tau) ---
            val = CdC(l1,l2) * conj(dB_az(:,l1)' * B(:,l2)) * AA(l1,l2);
            Omega(2*L+l1, 4*L+l2) = (2/sigma2) * real(val);
            Omega(4*L+l2, 2*L+l1) = (2/sigma2) * real(val);

            % --- Cross: (el_T, tau) ---
            val = CdC(l1,l2) * conj(dB_el(:,l1)' * B(:,l2)) * AA(l1,l2);
            Omega(3*L+l1, 4*L+l2) = (2/sigma2) * real(val);
            Omega(4*L+l2, 3*L+l1) = (2/sigma2) * real(val);

            % --- Cross: (az_R, az_T) ---
            val = cc * conj(dB_az(:,l2)' * B(:,l1)) * (dA_az(:,l1)' * A(:,l2));
            Omega(0*L+l1, 2*L+l2) = (2/sigma2) * real(val);
            Omega(2*L+l2, 0*L+l1) = (2/sigma2) * real(val);

            % --- Cross: (az_R, el_T) ---
            val = cc * conj(dB_el(:,l2)' * B(:,l1)) * (dA_az(:,l1)' * A(:,l2));
            Omega(0*L+l1, 3*L+l2) = (2/sigma2) * real(val);
            Omega(3*L+l2, 0*L+l1) = (2/sigma2) * real(val);

            % --- Cross: (el_R, az_T) ---
            val = cc * conj(dB_az(:,l2)' * B(:,l1)) * (dA_el(:,l1)' * A(:,l2));
            Omega(1*L+l1, 2*L+l2) = (2/sigma2) * real(val);
            Omega(2*L+l2, 1*L+l1) = (2/sigma2) * real(val);

            % --- Cross: (el_R, el_T) ---
            val = cc * conj(dB_el(:,l2)' * B(:,l1)) * (dA_el(:,l1)' * A(:,l2));
            Omega(1*L+l1, 3*L+l2) = (2/sigma2) * real(val);
            Omega(3*L+l2, 1*L+l1) = (2/sigma2) * real(val);
        end
    end

    % Scale matrix to improve conditioning before inversion
    d_vec = diag(Omega);
    d_vec(d_vec < eps) = 1; % Prevent divide by zero
    D_scale = diag(1 ./ sqrt(d_vec));
    
    Omega_scaled = D_scale * Omega * D_scale;
    
    disp('diag(Omega) values:');
    disp(diag(Omega)');
    
    % Regularize scaled matrix gently
    Omega_scaled = Omega_scaled + eye(5*L) * 1e-10;
    
    % Invert and descale
    CRB_scaled = inv(Omega_scaled);
    CRB_params = D_scale * CRB_scaled * D_scale;

    % Ensure CRB diagonal is positive
    for i = 1:5*L
        CRB_params(i,i) = max(0, real(CRB_params(i,i)));
    end

    % UT Position CRB via Jacobian (Eq. A32-A34)
    Jac_pR = zeros(3, 5*L);
    for l = 1:L
        pl = pl_all(:,l);
        dxy2  = (pl(1)-pR(1))^2 + (pl(2)-pR(2))^2 + eps;
        dRL   = dcR(l) + eps;
        denom_el = dRL * sqrt(max(eps, dRL^2-(pR(3)-pl(3))^2));

        % ∂θ^{az}_R/∂pR
        Jac_pR(1,(0)*L+l) = -(pl(2)-pR(2))/dxy2;
        Jac_pR(2,(0)*L+l) =  (pl(1)-pR(1))/dxy2;

        % ∂θ^{el}_R/∂pR
        Jac_pR(1,(1)*L+l) = -(pl(3)-pR(3))*(pl(1)-pR(1))/(denom_el + eps);
        Jac_pR(2,(1)*L+l) = -(pl(3)-pR(3))*(pl(2)-pR(2))/(denom_el + eps);
        Jac_pR(3,(1)*L+l) =  ((pl(1)-pR(1))^2+(pl(2)-pR(2))^2)/(denom_el + eps);

        % ∂τ/∂pR
        Jac_pR(1,(4)*L+l) = -(pl(1)-pR(1))/(c_speed*dRL);
        Jac_pR(2,(4)*L+l) = -(pl(2)-pR(2))/(c_speed*dRL);
        Jac_pR(3,(4)*L+l) = -(pl(3)-pR(3))/(c_speed*dRL);
    end

    FIM_UT = Jac_pR * Omega * Jac_pR' + eye(3)*1e-12;
    CRB_pR = pinv(FIM_UT);
    for i = 1:3, CRB_pR(i,i) = max(0, real(CRB_pR(i,i))); end

    % SP Position CRBs
    CRB_pl = zeros(3,3,L);
    for l = 1:L
        pl  = pl_all(:,l);
        dTL = dcT(l)+eps;
        dRL = dcR(l)+eps;
        dxy2_T = (pl(1)-pT(1))^2+(pl(2)-pT(2))^2+eps;
        denom_elT = dTL*sqrt(max(eps,dTL^2-(pT(3)-pl(3))^2));

        Jac_pl = zeros(3,5*L);

        % ∂θ^{az}_T/∂pl
        Jac_pl(1,(2)*L+l) =  (pl(2)-pT(2))/dxy2_T;
        Jac_pl(2,(2)*L+l) = -(pl(1)-pT(1))/dxy2_T;

        % ∂θ^{el}_T/∂pl
        Jac_pl(1,(3)*L+l) =  (pl(3)-pT(3))*(pl(1)-pT(1))/(denom_elT + eps);
        Jac_pl(2,(3)*L+l) =  (pl(3)-pT(3))*(pl(2)-pT(2))/(denom_elT + eps);
        Jac_pl(3,(3)*L+l) = -((pl(1)-pT(1))^2+(pl(2)-pT(2))^2)/(denom_elT + eps);

        % ∂τ/∂pl (both legs)
        Jac_pl(1,(4)*L+l) = (pl(1)-pT(1))/(c_speed*dTL) + (pl(1)-pR(1))/(c_speed*dRL);
        Jac_pl(2,(4)*L+l) = (pl(2)-pT(2))/(c_speed*dTL) + (pl(2)-pR(2))/(c_speed*dRL);
        Jac_pl(3,(4)*L+l) = (pl(3)-pT(3))/(c_speed*dTL) + (pl(3)-pR(3))/(c_speed*dRL);

        FIM_SP = Jac_pl * Omega * Jac_pl' + eye(3)*1e-12;
        CRB_pl(:,:,l) = pinv(FIM_SP);
        for i = 1:3, CRB_pl(i,i,l) = max(0, real(CRB_pl(i,i,l))); end
    end
end
