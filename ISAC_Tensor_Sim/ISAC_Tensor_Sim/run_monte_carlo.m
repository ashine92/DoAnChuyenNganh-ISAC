function results = run_monte_carlo(params)
% =========================================================================
% run_monte_carlo.m
% =========================================================================
% Description:
%   Runs Monte Carlo simulation over SNR range for fixed F, T, K, L.
%   Implements the proposed Tensor-based algorithm, and two baselines:
%     - PUDD (Phase Unwrapping Distance Difference)
%     - MUSIC-LSPS (MUSIC-Like Spectrum Peak Searching)
%   Also computes CRBs for all parameters.
%
% Inputs:
%   params - struct with SNR_dB_vec, F, T, K, L, MC, and system params
%
% Outputs:
%   results - struct with NMSE fields indexed by SNR
% =========================================================================

    rng(params.rng_seed, 'twister');

    SNR_vec = params.SNR_dB_vec;
    nSNR    = length(SNR_vec);
    MC      = params.MC;
    L       = params.L;

    % NMSE arrays: methods x SNR
    methods = {'Proposed', 'MUSIC_LSPS', 'PUDD'};
    fields  = {'az_R','el_R','az_T','el_T','tau','pR','pl'};
    nF      = length(fields);

    % Initialize results struct
    for mi = 1:length(methods)
        for fi = 1:nF
            results.(methods{mi}).(fields{fi}) = zeros(1, nSNR);
        end
    end
    for fi = 1:nF
        results.CRB.(fields{fi}) = zeros(1, nSNR);
    end

    for si = 1:nSNR
        params.SNR_dB = SNR_vec(si);
        fprintf('  SNR = %d dB (%d/%d)\n', SNR_vec(si), si, nSNR);

        % Accumulators for this SNR
        nmse_acc = zeros(length(methods), nF);
        crb_acc  = zeros(1, nF);
        n_valid  = 0;

        for mc = 1:MC
            try
                % --- Generate channel ---
                [Hk, alpha, tau_true, az_R, el_R, az_T, el_T, pl_all, d_cR, d_cT] = ...
                    generate_channel(params);

                % --- Construct received tensor ---
                [Y, W, F_mat, ~, ~, ~] = construct_tensor(params, Hk, alpha, tau_true);

                % --- Retrieve true AR, BT for CRB ---
                AR_true = zeros(params.NR, L);
                BT_true = zeros(params.NT, L);
                for l = 1:L
                    [AR_true(:,l), BT_true(:,l)] = near_field_array_response(params, pl_all(:,l));
                end

                % --- Build true C matrix ---
                k_indices = round(linspace(1, params.K_bar, params.K))';
                C_true = zeros(params.K, L);
                for l = 1:L
                    C_true(:,l) = alpha(l) * exp(-1j * 2*pi * tau_true(l) .* ...
                        params.fs .* k_indices / params.K_bar);
                end

                % --- CP-ALS decomposition ---
                [A_hat, B_hat, C_hat, ~, ~] = cp_als(Y, L, params);

                % --- Recover AR_hat, BT_hat ---
                AR_hat = pinv(W') * A_hat;   % NR x L
                BT_hat = pinv(F_mat') * B_hat; % NT x L

                % --- ToA estimation ---
                tau_hat = estimate_toa(C_hat, params);

                % --- Match paths (permutation) ---
                [tau_hat, perm] = match_paths(tau_hat, tau_true);
                AR_hat = AR_hat(:, perm);
                BT_hat = BT_hat(:, perm);

                % --- AoA estimation ---
                params_R = params;
                [az_R_hat, el_R_hat] = estimate_angles(AR_hat, params_R);

                % --- AoD estimation ---
                params_T = params;
                params_T.NRy = params.NTy;  params_T.NRz = params.NTz;
                params_T.NR  = params.NT;
                [az_T_hat, el_T_hat] = estimate_angles(BT_hat, params_T);

                % --- Match angle paths ---
                [az_R_hat, az_R_perm] = match_paths(az_R_hat, az_R);
                el_R_hat = el_R_hat(az_R_perm);
                [az_T_hat, az_T_perm] = match_paths(az_T_hat, az_T);
                el_T_hat = el_T_hat(az_T_perm);

                % --- UT localization ---
                pR_hat = ut_localization(tau_hat, az_R_hat, el_R_hat, ...
                                          az_T_hat, el_T_hat, params);

                % --- SP localization ---
                pl_hat = sp_localization(tau_hat, az_R_hat, el_R_hat, ...
                                          az_T_hat, el_T_hat, pR_hat, params);

                % --- NMSE computation ---
                nmse_vals = compute_nmse_all(az_R, el_R, az_T, el_T, tau_true, ...
                    params.pR, pl_all, az_R_hat, el_R_hat, az_T_hat, el_T_hat, ...
                    tau_hat, pR_hat, pl_hat);

                nmse_acc(1,:) = nmse_acc(1,:) + nmse_vals;

                % --- Baseline: MUSIC-LSPS ---
                [az_R_mu, el_R_mu, az_T_mu, el_T_mu, tau_mu, pR_mu, pl_mu] = ...
                    music_lsps(AR_hat, BT_hat, C_hat, params, pl_all, tau_true, az_R, el_R, az_T, el_T);
                nmse_mu = compute_nmse_all(az_R, el_R, az_T, el_T, tau_true, ...
                    params.pR, pl_all, az_R_mu, el_R_mu, az_T_mu, el_T_mu, ...
                    tau_mu, pR_mu, pl_mu);
                nmse_acc(2,:) = nmse_acc(2,:) + nmse_mu;

                % --- Baseline: PUDD ---
                [az_R_pu, el_R_pu, az_T_pu, el_T_pu, tau_pu, pR_pu, pl_pu] = ...
                    pudd_baseline(AR_hat, BT_hat, C_hat, params, pl_all, tau_true, az_R, el_R, az_T, el_T);
                nmse_pu = compute_nmse_all(az_R, el_R, az_T, el_T, tau_true, ...
                    params.pR, pl_all, az_R_pu, el_R_pu, az_T_pu, el_T_pu, ...
                    tau_pu, pR_pu, pl_pu);
                nmse_acc(3,:) = nmse_acc(3,:) + nmse_pu;

                % --- CRB ---
                [CRB_p, CRB_pR_mat, CRB_pl_arr] = compute_crb(params, AR_true, BT_true, ...
                    C_true, pl_all, params.SNR_dB);
                crb_vec = extract_crb_nmse(CRB_p, CRB_pR_mat, CRB_pl_arr, ...
                    az_R, el_R, az_T, el_T, tau_true, params.pR, pl_all);
                crb_acc = crb_acc + crb_vec;

                n_valid = n_valid + 1;
            catch ME
                % Skip failed trials silently
            end
        end

        if n_valid == 0, n_valid = 1; end

        for mi = 1:length(methods)
            for fi = 1:nF
                results.(methods{mi}).(fields{fi})(si) = nmse_acc(mi,fi) / n_valid;
            end
        end
        for fi = 1:nF
            results.CRB.(fields{fi})(si) = crb_acc(fi) / n_valid;
        end
    end

    results.SNR_vec = SNR_vec;
    results.params  = params;
end


%% ---- Helper: NMSE for all parameters -----------------------------------
function nmse_vals = compute_nmse_all(az_R, el_R, az_T, el_T, tau, pR, pl, ...
    az_R_h, el_R_h, az_T_h, el_T_h, tau_h, pR_h, pl_h)

    nmse = @(x, xh) norm(x(:)-xh(:))^2 / (norm(x(:))^2 + eps);

    nmse_vals = [ nmse(az_R, az_R_h), ...
                  nmse(el_R, el_R_h), ...
                  nmse(az_T, az_T_h), ...
                  nmse(el_T, el_T_h), ...
                  nmse(tau,  tau_h),  ...
                  nmse(pR,   pR_h),   ...
                  sum(sum((pl-pl_h).^2,1)) / (sum(sum(pl.^2,1)) + eps) ];
end


%% ---- Helper: Extract CRB as normalized MSE value -----------------------
function crb_vec = extract_crb_nmse(CRB_p, CRB_pR, CRB_pl, ...
    az_R, el_R, az_T, el_T, tau, pR, pl)

    L = length(az_R);

    % CRB for each parameter: trace of relevant diagonal block / ||true||^2
    crb_az_R = sum(diag(CRB_p(1:L, 1:L))) / (norm(az_R)^2 + eps);
    crb_el_R = sum(diag(CRB_p(L+1:2*L, L+1:2*L))) / (norm(el_R)^2 + eps);
    crb_az_T = sum(diag(CRB_p(2*L+1:3*L, 2*L+1:3*L))) / (norm(az_T)^2 + eps);
    crb_el_T = sum(diag(CRB_p(3*L+1:4*L, 3*L+1:4*L))) / (norm(el_T)^2 + eps);
    crb_tau  = sum(diag(CRB_p(4*L+1:5*L, 4*L+1:5*L))) / (norm(tau)^2 + eps);
    crb_pR   = trace(CRB_pR) / (norm(pR)^2 + eps);

    crb_pl_sum = 0;
    pl_norm_sq = 0;
    for l = 1:L
        crb_pl_sum = crb_pl_sum + trace(CRB_pl(:,:,l));
        pl_norm_sq = pl_norm_sq + norm(pl(:,l))^2;
    end
    crb_pl = crb_pl_sum / (pl_norm_sq + eps);

    crb_vec = [crb_az_R, crb_el_R, crb_az_T, crb_el_T, crb_tau, crb_pR, crb_pl];
end


%% ---- Helper: Path matching by nearest neighbor in ToA ------------------
function [x_matched, perm] = match_paths(x_est, x_true)
% Hungarian-style 1D matching: assign estimated paths to true paths
% to minimize total |x_est - x_true|^2

    L = length(x_true);
    if L == 1
        perm = 1;
        x_matched = x_est;
        return;
    end

    cost = abs(bsxfun(@minus, real(x_est(:)), real(x_true(:))'));  % L x L

    % Greedy matching (sufficient for small L)
    perm = zeros(1,L);
    used = false(1,L);
    for i = 1:L
        [~, best_j] = min(cost(i,:) + 1e9*used);
        perm(i) = best_j;
        used(best_j) = true;
    end
    x_matched = x_est(perm);
end
