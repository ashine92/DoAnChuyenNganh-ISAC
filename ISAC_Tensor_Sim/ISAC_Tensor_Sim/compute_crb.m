function [CRB_params, CRB_pR, CRB_pl] = compute_crb(params, AR, BT, C, pl_all, SNR_dB)
% =========================================================================
% compute_crb.m  —  Cramér-Rao Bound for channel parameters and positions
% =========================================================================
% Derives the FIM Ω(ϖ) ∈ R^{5L×5L} for
%   ϖ = [θ^{az}_R, θ^{el}_R, θ^{az}_T, θ^{el}_T, τ]
% and position CRBs via Jacobian transforms (Appendix B of paper).
%
% FIM block for same-type parameters, e.g. az_R (Eq. A19):
%   Ω_{az_R}(l1,l2) = (2/ς²) Re{ (dA_az_l1^H * dA_az_l2) * (CB_l1^H * CB_l2) }
%
% where dA_az_l = W^T * Λ^{az}_A * a_{R,l}  (Eq. A11)
% =========================================================================

    L      = params.L;
    NR     = params.NR;   NT = params.NT;
    NRy    = params.NRy;  NRz = params.NRz;
    NTy    = params.NTy;  NTz = params.NTz;
    K      = params.K;    K_bar = params.K_bar;
    F_size = params.F;    T_size = params.T;
    fs     = params.fs;   d = params.d;
    lambda = params.lambda; c = params.c;
    pT     = params.pT;   pR = params.pR;

    SNR_lin = 10^(SNR_dB/10);
    sigma2  = 1 / SNR_lin;

    % Subcarrier indices
    k_idx = round(linspace(1, K_bar, K))';   % K×1

    % DFT beamforming matrices
    W     = dftmtx(NR)/sqrt(NR);  W     = W(:,1:F_size);   % NR×F
    F_mat = dftmtx(NT)/sqrt(NT);  F_mat = F_mat(:,1:T_size); % NT×T

    % Effective factor matrices A=W'*AR, B=F'*BT
    A = W'  * AR;   % F×L
    B = F_mat' * BT; % T×L

    % KR products used in FIM
    CB = khatri_rao_crb(C, B);   % (K*T)×L
    CA = khatri_rao_crb(C, A);   % (K*F)×L
    BA = khatri_rao_crb(B, A);   % (T*F)×L

    % -------------------------------------------------------
    % Antenna index vectors for UT (NR×1 column vectors)
    % -------------------------------------------------------
    [NY_R, NZ_R] = meshgrid(-(NRy-1)/2:(NRy-1)/2, -(NRz-1)/2:(NRz-1)/2);
    ny_R = NY_R(:);   % NR×1
    nz_R = NZ_R(:);   % NR×1

    [NY_T, NZ_T] = meshgrid(-(NTy-1)/2:(NTy-1)/2, -(NTz-1)/2:(NTz-1)/2);
    ny_T = NY_T(:);   % NT×1
    nz_T = NZ_T(:);

    % -------------------------------------------------------
    % Recover true geometry from pl_all
    % -------------------------------------------------------
    az_R = zeros(L,1); el_R = zeros(L,1);
    az_T = zeros(L,1); el_T = zeros(L,1);
    tau  = zeros(L,1);
    dcR  = zeros(L,1); dcT = zeros(L,1);

    for l = 1:L
        pl      = pl_all(:,l);
        dcR(l)  = norm(pl - pR);
        dcT(l)  = norm(pl - pT);
        az_R(l) = atan2(pl(2)-pR(2), pl(1)-pR(1)) + pi;
        el_R(l) = acos(max(-1,min(1,(pR(3)-pl(3))/dcR(l))));
        az_T(l) = atan2(pl(2)-pT(2), pl(1)-pT(1));
        el_T(l) = acos(max(-1,min(1,(pT(3)-pl(3))/dcT(l))));
        tau(l)  = (dcR(l)+dcT(l))/c;
    end

    % -------------------------------------------------------
    % Build derivative matrices dA_az, dA_el, dB_az, dB_el, dC_tau
    % All are F×L or T×L or K×L  (NO nested functions — inline only)
    % -------------------------------------------------------
    dA_az  = zeros(F_size, L);
    dA_el  = zeros(F_size, L);
    dB_az  = zeros(T_size, L);
    dB_el  = zeros(T_size, L);
    dC_tau = zeros(K, L);

    for l = 1:L
        pl = pl_all(:,l);

        % Distances from SP l to each UT antenna (NR×1)
        ant_R  = pR + [zeros(1,NR); ny_R'*d; nz_R'*d];  % 3×NR
        dnR    = sqrt(sum((ant_R - pl).^2, 1))';          % NR×1

        % Distances from SP l to each BS antenna (NT×1)
        ant_T  = pT + [zeros(1,NT); ny_T'*d; nz_T'*d];  % 3×NT
        dnT    = sqrt(sum((ant_T - pl).^2, 1))';          % NT×1

        % ∂a_{R,l}/∂θ^{az}_{R,l}  (Eq. A11) — element-wise, NR×1 result
        % Λ^{az}_{nR} = j2πd*dcR*cos(az)*sin(el)/(λ*dnR) * ny_R  (diagonal elements)
        pfR_az = (1j*2*pi*d*dcR(l)*cos(az_R(l))*sin(el_R(l))/(lambda)) ...
                 ./ dnR;               % NR×1  (divide by dnR element-wise)
        daR_az = pfR_az .* ny_R .* AR(:,l);  % NR×1 element-wise
        dA_az(:,l) = W' * daR_az;            % F×1  ✓

        % ∂a_{R,l}/∂θ^{el}_{R,l}  (Eq. A14)
        pfR_el = (1j*2*pi*d*dcR(l)/lambda) ./ dnR;   % NR×1
        daR_el = pfR_el .* (sin(az_R(l))*cos(el_R(l))*ny_R ...
                            + sin(el_R(l))*nz_R) .* AR(:,l);
        dA_el(:,l) = W' * daR_el;            % F×1  ✓

        % ∂b_{T,l}/∂θ^{az}_{T,l}
        pfT_az = (1j*2*pi*d*dcT(l)*cos(az_T(l))*sin(el_T(l))/lambda) ...
                 ./ dnT;               % NT×1
        dbT_az = pfT_az .* ny_T .* BT(:,l);
        dB_az(:,l) = F_mat' * dbT_az;        % T×1  ✓

        % ∂b_{T,l}/∂θ^{el}_{T,l}
        pfT_el = (1j*2*pi*d*dcT(l)/lambda) ./ dnT;
        dbT_el = pfT_el .* (sin(az_T(l))*cos(el_T(l))*ny_T ...
                            + sin(el_T(l))*nz_T) .* BT(:,l);
        dB_el(:,l) = F_mat' * dbT_el;        % T×1  ✓

        % ∂c_l/∂τ_l  (Eq. A14) — element-wise on subcarrier index
        dC_tau(:,l) = -1j*2*pi*fs/K_bar * (k_idx-1) .* C(:,l);  % K×1
    end

    % -------------------------------------------------------
    % Build 5L×5L FIM (block diagonal, Eq. A19/A22)
    % Ω_{p,q}(l1,l2) = (2/ς²) Re{ (dM_p_l1^H * dM_p_l2) * (KR_p_l1^H * KR_p_l2) }
    % where M_p ∈ {dA_az, dA_el, dB_az, dB_el, dC_tau}
    %       KR_p ∈ {CB, CB, CA, CA, BA}  (matching mode unfolding per param)
    % -------------------------------------------------------
    D_mats = {dA_az, dA_el, dB_az, dB_el, dC_tau};  % 5 cells, each F/T/K × L
    KR_mats = {CB, CB, CA, CA, BA};                   % matching KR products

    Omega = zeros(5*L, 5*L);
    for p = 1:5
        Dp  = D_mats{p};   % F (or T or K) × L
        KRp = KR_mats{p};  % (K*T or K*F or T*F) × L
        for l1 = 1:L
            for l2 = 1:L
                row = (p-1)*L + l1;
                col = (p-1)*L + l2;
                Omega(row,col) = (2/sigma2) * real( ...
                    (Dp(:,l1)'*Dp(:,l2)) * (KRp(:,l1)'*KRp(:,l2)) );
            end
        end
    end

    % Add regularisation for invertibility
    Omega = Omega + eye(5*L)*1e-12;
    CRB_params = inv(Omega);

    % -------------------------------------------------------
    % UT Position CRB via Jacobian (Eq. A32-A34)
    % ∇_{pR} ϖ  is 3×5L
    % -------------------------------------------------------
    Jac_pR = zeros(3, 5*L);
    for l = 1:L
        pl = pl_all(:,l);
        dxy2  = (pl(1)-pR(1))^2 + (pl(2)-pR(2))^2 + eps;
        dRL   = dcR(l) + eps;
        denom_el = dRL * sqrt(max(eps, dRL^2-(pR(3)-pl(3))^2));

        % ∂θ^{az}_R/∂pR  (Eq. A34)
        Jac_pR(1,(0)*L+l) = (pR(2)-pl(2))/dxy2;
        Jac_pR(2,(0)*L+l) = (pl(1)-pR(1))/dxy2;
        Jac_pR(3,(0)*L+l) = 0;

        % ∂θ^{el}_R/∂pR
        Jac_pR(1,(1)*L+l) =  (pR(3)-pl(3))*(pR(1)-pl(1))/denom_el;
        Jac_pR(2,(1)*L+l) =  (pR(3)-pl(3))*(pR(2)-pl(2))/denom_el;
        Jac_pR(3,(1)*L+l) = -((pl(1)-pR(1))^2+(pl(2)-pR(2))^2)/denom_el;

        % ∂θ^{az}_T/∂pR = 0, ∂θ^{el}_T/∂pR = 0  (AoD ≠ f(pR))

        % ∂τ/∂pR
        Jac_pR(1,(4)*L+l) = (pR(1)-pl(1))/(c*dRL);
        Jac_pR(2,(4)*L+l) = (pR(2)-pl(2))/(c*dRL);
        Jac_pR(3,(4)*L+l) = (pR(3)-pl(3))/(c*dRL);
    end

    Omega_UT = Jac_pR * Omega * Jac_pR' + eye(3)*1e-12;
    CRB_pR   = inv(Omega_UT);

    % -------------------------------------------------------
    % SP Position CRBs
    % -------------------------------------------------------
    CRB_pl = zeros(3,3,L);
    for l = 1:L
        pl  = pl_all(:,l);
        dTL = dcT(l)+eps;
        dRL = dcR(l)+eps;
        dxy2_T = (pl(1)-pT(1))^2+(pl(2)-pT(2))^2+eps;
        denom_elT = dTL*sqrt(max(eps,dTL^2-(pT(3)-pl(3))^2));

        Jac_pl = zeros(3,5*L);

        % ∂θ^{az}_T/∂pl
        Jac_pl(1,(2)*L+l) = (pT(2)-pl(2))/dxy2_T;
        Jac_pl(2,(2)*L+l) = (pl(1)-pT(1))/dxy2_T;

        % ∂θ^{el}_T/∂pl
        Jac_pl(1,(3)*L+l) =  (pT(3)-pl(3))*(pl(1)-pT(1))/denom_elT;
        Jac_pl(2,(3)*L+l) =  (pT(3)-pl(3))*(pl(2)-pT(2))/denom_elT;
        Jac_pl(3,(3)*L+l) = -((pl(1)-pT(1))^2+(pl(2)-pT(2))^2)/denom_elT;

        % ∂τ/∂pl  (both legs)
        Jac_pl(1,(4)*L+l) = (pl(1)-pT(1))/(c*dTL) + (pl(1)-pR(1))/(c*dRL);
        Jac_pl(2,(4)*L+l) = (pl(2)-pT(2))/(c*dTL) + (pl(2)-pR(2))/(c*dRL);
        Jac_pl(3,(4)*L+l) = (pl(3)-pT(3))/(c*dTL) + (pl(3)-pR(3))/(c*dRL);

        Omega_SP = Jac_pl * Omega * Jac_pl' + eye(3)*1e-12;
        CRB_pl(:,:,l) = inv(Omega_SP);
    end
end


function KR = khatri_rao_crb(A, B)
% Khatri-Rao product: KR(:,l) = kron(A(:,l), B(:,l))
    [m, L] = size(A);  n = size(B,1);
    KR = reshape(bsxfun(@times, reshape(A,[m,1,L]), reshape(B,[1,n,L])), m*n, L);
end
