function results = run_monte_carlo(params)
% =========================================================================
% run_monte_carlo.m
% =========================================================================
% Runs Monte Carlo simulation over SNR range for fixed F, T, K, L.
% Implements: Proposed (Tensor-based), MUSIC-LSPS, PUDD, and CRB.
% =========================================================================

    rng(params.rng_seed, 'twister');

    SNR_vec = params.SNR_dB_vec;
    nSNR    = length(SNR_vec);
    MC      = params.MC;
    L       = params.L;

    methods = {'Proposed', 'MUSIC_LSPS', 'PUDD'};
    fields  = {'az_R','el_R','az_T','el_T','tau','pR','pl'};
    nF      = length(fields);

    % Initialize results
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

        nmse_all = nan(MC, length(methods), nF);
        norm_sq_all = nan(MC, nF);
        crb_all  = nan(MC, nF);
        crb_norm_sq_all = nan(MC, nF);
        n_valid  = 0;

        for mc = 1:MC
            try
                warning('off', 'MATLAB:rankDeficientMatrix');
                warning('off', 'MATLAB:nearlySingularMatrix');
                rng(mc); % Ensure same channel realization across SNRs for smooth curves
                % --- Generate channel ---
                [Hk, alpha, tau_true, az_R, el_R, az_T, el_T, pl_all, d_cR, d_cT] = ...
                    generate_channel(params);

                % --- Construct received tensor ---
                [Y, W, F_mat, ~, AR_true, BT_true] = ...
                    construct_tensor(params, Hk, alpha, tau_true, pl_all);

                % --- Build true C matrix ---
                k_indices = round(linspace(1, params.K_bar, params.K))';
                C_true = zeros(params.K, L);
                for l = 1:L
                    C_true(:,l) = alpha(l) * exp(-1j * 2*pi * tau_true(l) .* ...
                        params.fs .* k_indices / params.K_bar);
                end

                % --- CP-ALS decomposition ---
                [A_hat, B_hat, C_hat, ~, ~] = cp_als(Y, L, params);

                % --- FIX: Reorder factors to canonical order (by C-factor phase slopes) ---
                % CP decomposition has permutation invariance; reordering ensures consistent path assignment
                [A_hat, B_hat, C_hat, ~] = reorder_cp_factors(A_hat, B_hat, C_hat, params);

                % --- Recover AR_hat, BT_hat ---
                AR_hat = W * A_hat;    % NR x L
                BT_hat = F_mat * conj(B_hat); % NT x L

                % --- ToA estimation ---
                tau_hat = estimate_toa(C_hat, params);

                % --- AoA estimation ---
                params_R = params;
                [az_R_hat, el_R_hat] = estimate_angles(AR_hat, params_R);

                [perm] = match_paths_combined(az_R_hat, el_R_hat, tau_hat, ...
                                              az_R, el_R, tau_true);
                az_R_hat = az_R_hat(perm);
                el_R_hat = el_R_hat(perm);
                tau_hat  = tau_hat(perm);
                BT_hat   = BT_hat(:, perm);

                % --- AoD estimation ---
                params_T = params;
                params_T.NRy = params.NTy;
                params_T.NRz = params.NTz;
                params_T.NR  = params.NT;
                [az_T_hat, el_T_hat] = estimate_angles(BT_hat, params_T);

                % --- UT localization ---
                pR_hat = localize_ut(tau_hat, az_R_hat, el_R_hat, ...
                                     az_T_hat, el_T_hat, params);

                % --- SP localization ---
                pl_hat = localize_sps(tau_hat, az_R_hat, el_R_hat, ...
                                      az_T_hat, el_T_hat, pR_hat, params);

                % --- NMSE computation ---
                [nmse_vals, norms_sq] = compute_nmse(az_R, el_R, az_T, el_T, tau_true, ...
                    params.pR, pl_all, az_R_hat, el_R_hat, az_T_hat, el_T_hat, ...
                    tau_hat, pR_hat, pl_hat);

                nmse_all(mc, 1, :) = nmse_vals;
                norm_sq_all(mc, :) = norms_sq;

                % --- Baseline: MUSIC-LSPS ---
                if isfield(params, 'run_music') && params.run_music
                    [az_R_mu, el_R_mu, az_T_mu, el_T_mu, tau_mu, pR_mu, pl_mu] = ...
                        music_lsps(AR_hat, BT_hat, C_hat, params, pl_all, ...
                                   tau_true, az_R, el_R, az_T, el_T);
                    [nmse_mu, ~] = compute_nmse(az_R, el_R, az_T, el_T, tau_true, ...
                        params.pR, pl_all, az_R_mu, el_R_mu, az_T_mu, el_T_mu, ...
                        tau_mu, pR_mu, pl_mu);
                    if ~any(nmse_mu > 100)
                        nmse_all(mc, 2, :) = nmse_mu;
                    end
                end

                % --- Baseline: PUDD ---
                if isfield(params, 'run_pudd') && params.run_pudd
                    [az_R_pu, el_R_pu, az_T_pu, el_T_pu, tau_pu, pR_pu, pl_pu] = ...
                        pudd_baseline(AR_hat, BT_hat, C_hat, params, pl_all, ...
                                      tau_true, az_R, el_R, az_T, el_T);
                    [nmse_pu, ~] = compute_nmse(az_R, el_R, az_T, el_T, tau_true, ...
                        params.pR, pl_all, az_R_pu, el_R_pu, az_T_pu, el_T_pu, ...
                        tau_pu, pR_pu, pl_pu);
                    if ~any(nmse_pu > 100)
                        nmse_all(mc, 3, :) = nmse_pu;
                    end
                end

                % --- CRB ---
                if isfield(params, 'run_crb') && params.run_crb
                    [CRB_p, CRB_pR_mat, CRB_pl_arr] = compute_crb(params, ...
                        AR_true, BT_true, C_true, pl_all, params.SNR_dB);
                    [crb_vec, crb_norms_sq] = extract_crb_nmse(CRB_p, CRB_pR_mat, CRB_pl_arr, ...
                        az_R, el_R, az_T, el_T, tau_true, params.pR, pl_all);
                    crb_all(mc, :) = crb_vec;
                    crb_norm_sq_all(mc, :) = crb_norms_sq;
                end

                n_valid = n_valid + 1;

            catch ME
                % Print full error on first failure for debugging
                if mc <= 2
                    fprintf('    MC %d failed: %s\n', mc, ME.message);
                    for si2 = 1:length(ME.stack)
                        fprintf('      at %s (line %d)\n', ME.stack(si2).name, ME.stack(si2).line);
                    end
                end
            end
        end

        if n_valid == 0
            fprintf('    WARNING: no valid trials at SNR=%d dB\n', SNR_vec(si));
            n_valid = 1;
        else
            fprintf('    %d/%d valid trials\n', n_valid, MC);
        end

        for mi = 1:length(methods)
            for fi = 1:nF
                % Median MSE / Mean Squared Norm * 100
                med_mse = nanmedian(nmse_all(:, mi, fi));
                mean_norm_sq = nanmean(norm_sq_all(:, fi));
                results.(methods{mi}).(fields{fi})(si) = (med_mse / mean_norm_sq) * 100;
            end
        end
        for fi = 1:nF
            med_crb = nanmedian(crb_all(:, fi));
            mean_crb_norm_sq = nanmean(crb_norm_sq_all(:, fi));
            results.CRB.(fields{fi})(si) = (med_crb / mean_crb_norm_sq) * 100;
        end
    end

    results.SNR_vec = SNR_vec;
    results.params  = params;
end


%% ---- Helper: Extract CRB as normalized MSE value -----------------------
function [crb_vec, norms_sq] = extract_crb_nmse(CRB_p, CRB_pR, CRB_pl, ...
    az_R, el_R, az_T, el_T, tau, pR, pl)

    L = length(az_R);
    dim = size(CRB_p, 1);

    % Guard against dimension mismatch
    if dim < 5*L
        crb_vec = zeros(1, 7);
        norms_sq = zeros(1, 7);
        return;
    end

    crb_az_R = sum(diag(CRB_p(1:L, 1:L)));
    crb_el_R = sum(diag(CRB_p(L+1:2*L, L+1:2*L)));
    crb_az_T = sum(diag(CRB_p(2*L+1:3*L, 2*L+1:3*L)));
    crb_el_T = sum(diag(CRB_p(3*L+1:4*L, 3*L+1:4*L)));
    crb_tau  = sum(diag(CRB_p(4*L+1:5*L, 4*L+1:5*L)));
    crb_pR   = trace(CRB_pR);

    crb_pl_sum = 0;
    pl_norm_sq = 0;
    for l = 1:L
        crb_pl_sum = crb_pl_sum + trace(CRB_pl(:,:,l));
        pl_norm_sq = pl_norm_sq + norm(pl(:,l))^2;
    end
    crb_pl = crb_pl_sum;

    crb_vec = [crb_az_R, crb_el_R, crb_az_T, crb_el_T, crb_tau, crb_pR, crb_pl];
    norms_sq = [norm(az_R)^2, norm(el_R)^2, norm(az_T)^2, norm(el_T)^2, norm(tau)^2, norm(pR)^2, pl_norm_sq];
end


%% ---- Helper: Combined path matching ------------------------------------
function perm = match_paths_combined(az_est, el_est, tau_est, az_true, el_true, tau_true)
% Match estimated paths to true paths using combined cost on angles and ToA.
% Uses Hungarian-style greedy matching.

    L = length(az_true);
    if L == 1
        perm = 1;
        return;
    end

    % Normalized cost matrix (combine angle and ToA errors)
    cost_az  = abs(bsxfun(@minus, az_est(:), az_true(:)'));
    cost_el  = abs(bsxfun(@minus, el_est(:), el_true(:)'));
    cost_tau = abs(bsxfun(@minus, tau_est(:), tau_true(:)'));

    % Normalize each cost by its range
    range_az  = max(eps, max(cost_az(:)));
    range_el  = max(eps, max(cost_el(:)));
    range_tau = max(eps, max(cost_tau(:)));

    cost = cost_az/range_az + cost_el/range_el + cost_tau/range_tau;

    % Greedy matching (true j -> est i)
    perm = zeros(1, L);
    used = false(1, L);
    for j = 1:L
        c = cost(:,j);
        c(used) = inf;
        [~, best_i] = min(c);
        perm(j) = best_i;
        used(best_i) = true;
    end
end
