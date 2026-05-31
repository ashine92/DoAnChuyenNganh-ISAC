function [az_R_hat, el_R_hat, az_T_hat, el_T_hat, tau_hat, pR_hat, pl_hat] = ...
    pudd_baseline(AR_hat, BT_hat, C_hat, params, pl_all, tau_true, az_R, el_R, az_T, el_T)
% =========================================================================
% pudd_baseline.m — Phase Unwrapping Distance Difference (PUDD)
% =========================================================================
% Reference: [39] Podkurkov et al., IEEE TSP 2021
%
% Key difference from proposed method:
%   - Uses same ALS factor matrices but estimates distances via phase
%     unwrapping across antenna pairs
%   - Distance estimation accumulates errors, degrading localization
%   - Angle estimation same as proposed (ESPRIT-like from CP factors)
%   - Localization uses distances rather than purely geometric approach
% =========================================================================

    L = params.L;
    d = params.d;
    lambda = params.lambda;
    c_speed = params.c;
    NRy = params.NRy; NRz = params.NRz;
    NTy = params.NTy; NTz = params.NTz;

    % --- Same angle estimates as proposed ---
    params_R = params;
    [az_R_hat, el_R_hat] = estimate_angles(AR_hat, params_R);

    params_T = params;
    params_T.NRy = NTy; params_T.NRz = NTz; params_T.NR = params.NT;
    [az_T_hat, el_T_hat] = estimate_angles(BT_hat, params_T);

    % --- ToA from C_hat ---
    tau_hat = estimate_toa(C_hat, params);

    % --- Match paths ---
    cost = abs(bsxfun(@minus, az_R_hat(:), az_R(:)')) + ...
           abs(bsxfun(@minus, el_R_hat(:), el_R(:)'));
    perm = greedy_match_pudd(cost);
    az_R_hat = az_R_hat(perm);
    el_R_hat = el_R_hat(perm);
    az_T_hat = az_T_hat(perm);
    el_T_hat = el_T_hat(perm);
    tau_hat  = tau_hat(perm);

    % --- PUDD distance estimation (introduces accumulation errors) ---
    % Estimate distance from phase curvature across antenna pairs
    % The near-field quadratic phase term gives: Φ_quadratic ∝ d²/(2*λ*distance)
    % Phase unwrapping error accumulates over multiple antenna pairs
    d_R_est = zeros(L,1);
    d_T_est = zeros(L,1);

    AR_hat_perm = AR_hat(:, perm);
    BT_hat_perm = BT_hat(:, perm);

    for l = 1:L
        % Estimate distance from quadratic phase (PUDD method)
        aR_l = AR_hat_perm(:,l);
        aR_mat = reshape(aR_l, NRz, NRy);

        % Phase difference along y-axis
        diff_y = aR_mat(:, 2:end) .* conj(aR_mat(:, 1:end-1));
        phase_y = angle(diff_y);

        % Second-order phase difference (curvature) gives distance info
        if NRy > 2
            diff2_y = phase_y(:, 2:end) - phase_y(:, 1:end-1);
            curvature_y = mean(abs(diff2_y(:)));
            if curvature_y > eps
                d_R_est(l) = (2*pi*d^2) / (lambda * curvature_y);
            else
                d_R_est(l) = c_speed * tau_hat(l) / 2;
            end
        else
            d_R_est(l) = c_speed * tau_hat(l) / 2;
        end

        % BS side
        bT_l = BT_hat_perm(:,l);
        bT_mat = reshape(bT_l, NTz, NTy);
        diff_y_T = bT_mat(:, 2:end) .* conj(bT_mat(:, 1:end-1));
        phase_y_T = angle(diff_y_T);

        if NTy > 2
            diff2_y_T = phase_y_T(:, 2:end) - phase_y_T(:, 1:end-1);
            curv_T = mean(abs(diff2_y_T(:)));
            if curv_T > eps
                d_T_est(l) = (2*pi*d^2) / (lambda * curv_T);
            else
                d_T_est(l) = c_speed * tau_hat(l) / 2;
            end
        else
            d_T_est(l) = c_speed * tau_hat(l) / 2;
        end
    end

    % PUDD localization uses distances + angles
    % SP = pT + d_T * gT  and  SP = pR + d_R * gR
    % This accumulates errors from both distance and angle estimation
    pT = params.pT;
    pR_hat = zeros(3, 1);
    weight_sum = 0;

    for l = 1:L
        az_T_l = az_T_hat(l); el_T_l = el_T_hat(l);
        az_R_l = az_R_hat(l); el_R_l = el_R_hat(l);

        gT = [cos(az_T_l)*sin(el_T_l); sin(az_T_l)*sin(el_T_l); cos(el_T_l)];
        gR = [cos(az_R_l)*sin(el_R_l); sin(az_R_l)*sin(el_R_l); cos(el_R_l)];
        gT = gT / (norm(gT)+eps);
        gR = gR / (norm(gR)+eps);

        % SP position from BS side
        pl_est_T = pT + d_T_est(l) * gT;
        % UT = SP - d_R * gR
        pR_est_l = pl_est_T - d_R_est(l) * gR;

        w = 1 / (d_R_est(l) + d_T_est(l) + eps);  % weight by inverse distance
        pR_hat = pR_hat + w * pR_est_l;
        weight_sum = weight_sum + w;
    end
    pR_hat = pR_hat / (weight_sum + eps);

    % SP localization
    pl_hat = zeros(3, L);
    for l = 1:L
        gT = [cos(az_T_hat(l))*sin(el_T_hat(l)); sin(az_T_hat(l))*sin(el_T_hat(l)); cos(el_T_hat(l))];
        gT = gT / (norm(gT)+eps);
        pl_hat(:,l) = pT + d_T_est(l) * gT;
    end
end


function perm = greedy_match_pudd(cost)
    L = size(cost, 1);
    perm = zeros(1, L);
    used = false(1, L);
    for j = 1:L
        c = cost(:,j); c(used) = inf;
        [~, i] = min(c); perm(j) = i; used(i) = true;
    end
end
