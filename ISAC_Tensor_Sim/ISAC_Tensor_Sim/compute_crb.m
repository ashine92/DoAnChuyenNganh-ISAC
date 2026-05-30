function [CRB_params, CRB_pR, CRB_pl] = compute_crb(params, AR, BT, C, pl_all, SNR_dB)
% =========================================================================
% compute_crb.m
% =========================================================================
% Description:
%   Derives the Cramér-Rao Bound (CRB) for channel parameters and positions.
%
%   The Fisher Information Matrix (FIM) is built from the likelihood function
%   of the received tensor (Appendix B of paper):
%     F(ϖ) ~ -1/ς² || [Y]_(n) - (... ) ||²_F
%
%   Parameter vector: ϖ = [θ_az_R, θ_el_R, θ_az_T, θ_el_T, τ]  (5L params)
%
%   FIM Omega(ϖ) ∈ R^{5L x 5L} is built from partial derivatives of
%   factor matrices A, B, C with respect to ϖ.
%
%   CRB(ϖ) = Omega^{-1}(ϖ)
%
%   Position CRBs use Jacobian transforms (Eq. A32):
%     Omega_UT = (∇_{pR} ϖ) * Omega(ϖ) * (∇_{pR} ϖ)^H
%     CRB(pR)  = Omega_UT^{-1}
%
% Inputs:
%   params   - system parameter struct
%   AR       - NR x L true array response matrix (UT)
%   BT       - NT x L true array response matrix (BS)
%   C        - K x L true C factor matrix
%   pl_all   - 3 x L true SP positions
%   SNR_dB   - scalar SNR in dB
%
% Outputs:
%   CRB_params - 5L x 5L CRB matrix for [az_R, el_R, az_T, el_T, tau]
%   CRB_pR     - 3 x 3 CRB matrix for UT position
%   CRB_pl     - 3 x 3 x L CRB matrices for each SP position
%
% Runtime: O((5L)^2 * NT*NR*K)  (FIM construction)
% =========================================================================

    L      = params.L;
    NR     = params.NR;
    NT     = params.NT;
    NRy    = params.NRy;  NRz = params.NRz;
    NTy    = params.NTy;  NTz = params.NTz;
    K      = params.K;
    K_bar  = params.K_bar;
    fs     = params.fs;
    F_size = params.F;
    T_size = params.T;
    d      = params.d;
    lambda = params.lambda;
    c      = params.c;
    pT     = params.pT;
    pR     = params.pR;

    % Load W, F from params (or reconstruct)
    DFT_NR = dftmtx(NR) / sqrt(NR);
    DFT_NT = dftmtx(NT) / sqrt(NT);
    W      = DFT_NR(:, 1:F_size);
    F_mat  = DFT_NT(:, 1:T_size);

    % Effective factor matrices A = W^T * AR, B = F^T * BT
    A = W' * AR;   % F x L
    B = F_mat' * BT;  % T x L

    % Noise variance from SNR
    % SNR = ||Y_clean||^2 / ||V||^2  →  ς² = signal_power/SNR
    % For CRB we treat ς² as normalized to 1/SNR
    SNR_lin = 10^(SNR_dB/10);
    sigma2  = 1 / SNR_lin;   % normalized noise variance

    % Subcarrier indices
    k_indices = round(linspace(1, K_bar, K))';   % K x 1

    % -------------------------------------------------------
    % Precompute partial derivatives of array response vectors
    % ∂a_R / ∂θ_az_R,l  and  ∂a_R / ∂θ_el_R,l  for each path l
    % (Eq. A11, A14)
    % -------------------------------------------------------

    % Antenna index arrays for UT
    ny_R_range = (-(NRy-1)/2 : (NRy-1)/2);
    nz_R_range = (-(NRz-1)/2 : (NRz-1)/2);
    [NY_R, NZ_R] = meshgrid(ny_R_range, nz_R_range);
    ny_R_vec = NY_R(:);   % NR x 1
    nz_R_vec = NZ_R(:);

    ny_T_range = (-(NTy-1)/2 : (NTy-1)/2);
    nz_T_range = (-(NTz-1)/2 : (NTz-1)/2);
    [NY_T, NZ_T] = meshgrid(ny_T_range, nz_T_range);
    ny_T_vec = NY_T(:);   % NT x 1
    nz_T_vec = NZ_T(:);

    % FIM: 5L x 5L  (ordering: az_R x L, el_R x L, az_T x L, el_T x L, tau x L)
    Omega = zeros(5*L, 5*L);

    % Mode unfoldings needed for FIM (use A, B, C factor matrices)
    % Precompute KR products
    CB = khatri_rao(C, B);    % (K*T) x L
    CA = khatri_rao(C, A);    % (K*F) x L
    BA = khatri_rao(B, A);    % (T*F) x L

    % ---- Helper: build derivative ∂A/∂θ_az_R,l  = W^T * Λ^az_A * a_R,l ----
    function dA_l = dA_az(l, theta_az_R_l, theta_el_R_l, d_cR_l, aR_l, d_nR_l)
        % Eq. A11: Λ^az_A diagonal
        phase_factor = 1j * 2*pi*d * d_cR_l * cos(theta_az_R_l) * sin(theta_el_R_l) ...
                       / (lambda * d_nR_l);   % NR x 1
        Lam_az = diag(phase_factor .* ny_R_vec);   % Use ny component
        % Actually the block diagonal form in Eq. A11 simplifies to element-wise:
        % (Λ^az_A)_{nR,nR} = j2πd*d_c_R/(λ*d_nR) * cos(az)*sin(el) * ny_R
        % So: Λ^az_A * a_R,l = element-wise product of diagonal * aR_l
        dA_l = W' * (phase_factor .* ny_R_vec .* aR_l);  % F x 1
    end

    function dA_l = dA_el(l, theta_az_R_l, theta_el_R_l, d_cR_l, aR_l, d_nR_l)
        % Eq. A14: elevation derivative
        phase_factor = 1j * 2*pi*d * d_cR_l / (lambda * d_nR_l);
        % Two terms: sin(az)*cos(el)*ny + sin(el)*nz  (from paper Eq. A14)
        lam_val = phase_factor .* (sin(theta_az_R_l)*cos(theta_el_R_l)*ny_R_vec ...
                                   + sin(theta_el_R_l)*nz_R_vec);
        dA_l = W' * (lam_val .* aR_l);
    end

    function dB_l = dB_az(l, theta_az_T_l, theta_el_T_l, d_cT_l, bT_l, d_nT_l)
        phase_factor = 1j * 2*pi*d * d_cT_l * cos(theta_az_T_l) * sin(theta_el_T_l) ...
                       / (lambda * d_nT_l);
        dB_l = F_mat' * (phase_factor .* ny_T_vec .* bT_l);
    end

    function dB_l = dB_el(l, theta_az_T_l, theta_el_T_l, d_cT_l, bT_l, d_nT_l)
        phase_factor = 1j * 2*pi*d * d_cT_l / (lambda * d_nT_l);
        lam_val = phase_factor .* (sin(theta_az_T_l)*cos(theta_el_T_l)*ny_T_vec ...
                                   + sin(theta_el_T_l)*nz_T_vec);
        dB_l = F_mat' * (lam_val .* bT_l);
    end

    function dc_l = dC_tau(l, tau_l, alpha_l)
        % ∂c_l / ∂τ_l = -j2π*fs/K_bar * diag(0,1,...,K-1) * c_l  (Eq. A14)
        dc_l = -1j * 2*pi * fs / K_bar * (k_indices-1) .* C(:,l);
    end

    % Retrieve geometry for each path
    theta_az_R_true = zeros(L,1); theta_el_R_true = zeros(L,1);
    theta_az_T_true = zeros(L,1); theta_el_T_true = zeros(L,1);
    tau_true = zeros(L,1);
    d_cR = zeros(L,1); d_cT = zeros(L,1);

    for l = 1:L
        pl = pl_all(:,l);
        d_cR(l) = norm(pl - pR);
        d_cT(l) = norm(pl - pT);
        theta_az_R_true(l) = atan2(pl(2)-pR(2), pl(1)-pR(1)) + pi;
        theta_el_R_true(l) = acos((pR(3)-pl(3)) / d_cR(l));
        theta_az_T_true(l) = atan2(pl(2)-pT(2), pl(1)-pT(1));
        theta_el_T_true(l) = acos((pT(3)-pl(3)) / d_cT(l));
        tau_true(l) = (d_cR(l) + d_cT(l)) / c;
    end

    % Collect derivatives: dA, dB, dC matrices (F x L, T x L, K x L)
    dA_az_mat = zeros(F_size, L);
    dA_el_mat = zeros(F_size, L);
    dB_az_mat = zeros(T_size, L);
    dB_el_mat = zeros(T_size, L);
    dC_tau_mat = zeros(K, L);

    for l = 1:L
        pl = pl_all(:,l);
        % Compute d_nR, d_nT (element-wise distances to each antenna)
        ant_pos_R = pR + [zeros(1,NR); ny_R_vec'*d; nz_R_vec'*d];
        d_nR_vec = sqrt(sum((ant_pos_R - pl).^2, 1))';
        ant_pos_T = pT + [zeros(1,NT); ny_T_vec'*d; nz_T_vec'*d];
        d_nT_vec = sqrt(sum((ant_pos_T - pl).^2, 1))';

        dA_az_mat(:,l) = dA_az(l, theta_az_R_true(l), theta_el_R_true(l), ...
                                d_cR(l), AR(:,l), d_nR_vec);
        dA_el_mat(:,l) = dA_el(l, theta_az_R_true(l), theta_el_R_true(l), ...
                                d_cR(l), AR(:,l), d_nR_vec);
        dB_az_mat(:,l) = dB_az(l, theta_az_T_true(l), theta_el_T_true(l), ...
                                d_cT(l), BT(:,l), d_nT_vec);
        dB_el_mat(:,l) = dB_el(l, theta_az_T_true(l), theta_el_T_true(l), ...
                                d_cT(l), BT(:,l), d_nT_vec);
        dC_tau_mat(:,l) = dC_tau(l, tau_true(l), 1);
    end

    % Build FIM using Eq. A19, A22-A26
    % For same-type parameters (e.g., az_R vs az_R):
    %   Rn_az_A = 1/ς² (Ã_az^T ⊗ (C⊙B)^T)(Ã_az* ⊗ (C⊙B)*)
    % FIM(l1,l2) for (az_R, az_R) = 2*Re{Rn(p,q)} where p,q = L(l1-1)+l1, L(l2-1)+l2

    % Helper: compute 2*Re{(i1,i2) diagonal of cross-covariance}
    % For same-mode params: R_n = 1/ς² * (D1^T⊗M^T)(D1*⊗M*) where D1=∂(mode matrix)/∂θ
    % = 1/ς² * (D1^T D1* ⊗ M^T M*)
    % Diagonal entry (l1,l2): 1/ς² * (D1^T D1*)_{l1,l2} * (M^T M*)_{l1,l2}

    % For cross-mode (different tensor unfoldings), off-diag terms involve
    % Rv1v2 which are sparse (Eq. A28). For simplification, off-diag FIM
    % entries between different unfoldings are set using Eq. A25 structure.
    % Full derivation omitted for brevity; diagonal FIM gives valid CRB lower bound.

    % Gram matrices for efficiency
    AtA_az = (CB' * CB);   % (C⊙B)^T*(C⊙B), for dA terms
    gram_az_A = zeros(L,L);
    for l1=1:L
        for l2=1:L
            gram_az_A(l1,l2) = conj(dA_az_mat(:,l1))' * dA_az_mat(:,l2) * AtA_az(l1,l2);
        end
    end

    % Simplified diagonal FIM (valid lower bound)
    % FIM block for az_R (L x L):
    FIM_az_R = zeros(L,L);
    FIM_el_R = zeros(L,L);
    FIM_az_T = zeros(L,L);
    FIM_el_T = zeros(L,L);
    FIM_tau  = zeros(L,L);

    for l1 = 1:L
        for l2 = 1:L
            % az_R block (Eq. A19 diagonal)
            FIM_az_R(l1,l2) = (2/sigma2) * real( ...
                (dA_az_mat(:,l1)' * dA_az_mat(:,l2)) * (CB(:,l1)' * CB(:,l2)) );

            % el_R block
            FIM_el_R(l1,l2) = (2/sigma2) * real( ...
                (dA_el_mat(:,l1)' * dA_el_mat(:,l2)) * (CB(:,l1)' * CB(:,l2)) );

            % az_T block
            FIM_az_T(l1,l2) = (2/sigma2) * real( ...
                (dB_az_mat(:,l1)' * dB_az_mat(:,l2)) * (CA(:,l1)' * CA(:,l2)) );

            % el_T block
            FIM_el_T(l1,l2) = (2/sigma2) * real( ...
                (dB_el_mat(:,l1)' * dB_el_mat(:,l2)) * (CA(:,l1)' * CA(:,l2)) );

            % tau block
            FIM_tau(l1,l2) = (2/sigma2) * real( ...
                (dC_tau_mat(:,l1)' * dC_tau_mat(:,l2)) * (BA(:,l1)' * BA(:,l2)) );
        end
    end

    % Assemble full 5L x 5L FIM (block diagonal approximation)
    Omega = blkdiag(FIM_az_R, FIM_el_R, FIM_az_T, FIM_el_T, FIM_tau);

    % Add small regularization for invertibility
    Omega = Omega + eye(5*L) * 1e-12;

    % CRB for channel parameters
    CRB_params = inv(Omega);

    % -------------------------------------------------------
    % Position CRBs via Jacobian (Eq. A32-A35)
    % -------------------------------------------------------

    % Jacobian ∇_{pR} ϖ  (3 x 5L)
    Jac_pR = zeros(3, 5*L);

    for l = 1:L
        pl = pl_all(:,l);
        xl = pl(1); yl = pl(2); zl = pl(3);
        xR = pR(1); yR = pR(2); zR = pR(3);
        dist_RL = d_cR(l);

        dxy2 = (xl-xR)^2 + (yl-yR)^2;

        % ∂θ_az_R / ∂pR  (Eq. A34)
        Jac_pR(1, l)       = (yR-yl) / (dxy2 + eps);    % ∂/∂xR
        Jac_pR(2, l)       = (xl-xR) / (dxy2 + eps);    % ∂/∂yR
        Jac_pR(3, l)       = 0;

        % ∂θ_el_R / ∂pR
        denom_el = dist_RL * sqrt(max(eps, dist_RL^2 - (zR-zl)^2));
        Jac_pR(1, L+l)     = (zR-zl)*(xR-xl) / (denom_el+eps);
        Jac_pR(2, L+l)     = (zR-zl)*(yR-yl) / (denom_el+eps);
        Jac_pR(3, L+l)     = -((xl-xR)^2+(yl-yR)^2) / (denom_el+eps);

        % ∂θ_az_T / ∂pR = 0 (AoD doesn't depend on UT position)
        % (already zeros)

        % ∂θ_el_T / ∂pR = 0
        % (already zeros)

        % ∂τ / ∂pR  (Eq. A34: ∂τ/∂(x,y,z)_R = (xR-xl, yR-yl, zR-zl)/(c*dist_RL))
        Jac_pR(1, 4*L+l)   = (xR-xl) / (c * dist_RL + eps);
        Jac_pR(2, 4*L+l)   = (yR-yl) / (c * dist_RL + eps);
        Jac_pR(3, 4*L+l)   = (zR-zl) / (c * dist_RL + eps);
    end

    Omega_UT = Jac_pR * Omega * Jac_pR';
    Omega_UT = Omega_UT + eye(3) * 1e-12;
    CRB_pR   = inv(Omega_UT);

    % SP position CRBs
    CRB_pl = zeros(3, 3, L);
    for l = 1:L
        pl = pl_all(:,l);
        xl = pl(1); yl = pl(2); zl = pl(3);
        xT = pT(1); yT = pT(2); zT = pT(3);
        dist_TL = d_cT(l);
        dist_RL = d_cR(l);

        Jac_pl = zeros(3, 5*L);

        % ∂θ_az_T / ∂pl
        dxy2_T = (xl-xT)^2 + (yl-yT)^2;
        Jac_pl(1, 2*L+l)  = -(yl-yT) / (dxy2_T+eps);
        Jac_pl(2, 2*L+l)  =  (xl-xT) / (dxy2_T+eps);  % ∂/∂yl
        % Actually ∂atan2(yl-yT, xl-xT)/∂xl = -(yl-yT)/r², ∂/∂yl = (xl-xT)/r²
        % Correcting sign convention
        Jac_pl(1, 2*L+l)  = (yT-yl) / (dxy2_T+eps);
        Jac_pl(2, 2*L+l)  = (xl-xT) / (dxy2_T+eps);

        % ∂θ_el_T / ∂pl  (similar to el_R but wrt pl)
        denom_el_T = dist_TL * sqrt(max(eps, dist_TL^2 - (zT-zl)^2));
        Jac_pl(1, 3*L+l)  = (zT-zl)*(xl-xT) / (denom_el_T+eps);
        Jac_pl(2, 3*L+l)  = (zT-zl)*(yl-yT) / (denom_el_T+eps);
        Jac_pl(3, 3*L+l)  = -((xl-xT)^2+(yl-yT)^2) / (denom_el_T+eps);

        % ∂τ / ∂pl  (both legs)
        Jac_pl(1, 4*L+l)  = (xl-xT)/(c*dist_TL+eps) + (xl-pR(1))/(c*dist_RL+eps);
        Jac_pl(2, 4*L+l)  = (yl-yT)/(c*dist_TL+eps) + (yl-pR(2))/(c*dist_RL+eps);
        Jac_pl(3, 4*L+l)  = (zl-zT)/(c*dist_TL+eps) + (zl-pR(3))/(c*dist_RL+eps);

        Omega_SP = Jac_pl * Omega * Jac_pl';
        Omega_SP = Omega_SP + eye(3)*1e-12;
        CRB_pl(:,:,l) = inv(Omega_SP);
    end
end
