function results = run_monte_carlo_vs_K(params)
% =========================================================================
% run_monte_carlo_vs_K.m
% =========================================================================
% Runs Monte Carlo simulation over K (subcarrier count) range
% for fixed SNR, F, T, L. Used for Figure 5.
% =========================================================================

    rng(params.rng_seed, 'twister');

    N_vec = params.N_vec;
    nN    = length(N_vec);
    MC    = params.MC;
    L     = params.L;

    methods = {'Proposed', 'MUSIC_LSPS', 'PUDD'};
    fields  = {'az_R','el_R','az_T','el_T','tau','pR','pl'};

    for mi = 1:length(methods)
        for fi = 1:length(fields)
            results.(methods{mi}).(fields{fi}) = zeros(1, nN);
        end
    end
    for fi = 1:length(fields)
        results.CRB.(fields{fi}) = zeros(1, nN);
    end

    for ni = 1:nN
        N1D = N_vec(ni);
        params.NTy = N1D; params.NTz = N1D; params.NT = N1D^2;
        params.NRy = N1D; params.NRz = N1D; params.NR = N1D^2;
        params.F = N1D^2;
        params.T = N1D^2;
        
        fprintf('  N = %dx%d=%d (%d/%d)\n', N1D, N1D, N1D^2, ni, nN);

        nmse_all = nan(MC, length(methods), length(fields));
        norm_sq_all = nan(MC, length(fields));
        crb_all  = nan(MC, length(fields));
        crb_norm_sq_all = nan(MC, length(fields));
        n_valid  = 0;

        for mc = 1:MC
            try
                warning('off', 'MATLAB:rankDeficientMatrix');
                warning('off', 'MATLAB:nearlySingularMatrix');
                rng(mc);
                [Hk, alpha, tau_true, az_R, el_R, az_T, el_T, pl_all, ~, ~] = ...
                    generate_channel(params);
                [Y, W, F_mat, ~, AR_true, BT_true] = ...
                    construct_tensor(params, Hk, alpha, tau_true, pl_all);

                k_indices = round(linspace(1, params.K_bar, params.K))';
                C_true = zeros(params.K, L);
                for l = 1:L
                    C_true(:,l) = alpha(l) * exp(-1j*2*pi*tau_true(l).*params.fs.*k_indices/params.K_bar);
                end

                [A_hat, B_hat, C_hat, ~, ~] = cp_als(Y, L, params);

                % Recover AR, BT
                AR_hat = W * A_hat;
                BT_hat = F_mat * conj(B_hat);

                % ToA estimation
                tau_hat = estimate_toa(C_hat, params);

                % AoA estimation
                params_R = params;
                [az_R_hat, el_R_hat] = estimate_angles(AR_hat, params_R);

                % Match paths
                perm = match_paths_combined(az_R_hat, el_R_hat, tau_hat, ...
                                            az_R, el_R, tau_true);
                az_R_hat = az_R_hat(perm);
                el_R_hat = el_R_hat(perm);
                tau_hat  = tau_hat(perm);
                BT_hat   = BT_hat(:, perm);

                % AoD estimation
                params_T = params;
                params_T.NRy = params.NTy;
                params_T.NRz = params.NTz;
                params_T.NR  = params.NT;
                [az_T_hat, el_T_hat] = estimate_angles(BT_hat, params_T);

                % Localization
                pR_hat = localize_ut(tau_hat, az_R_hat, el_R_hat, ...
                                     az_T_hat, el_T_hat, params);
                pl_hat = localize_sps(tau_hat, az_R_hat, el_R_hat, ...
                                      az_T_hat, el_T_hat, pR_hat, params);

                % NMSE - Proposed
                [nmse_vals, norms_sq] = compute_nmse(az_R, el_R, az_T, el_T, tau_true, ...
                    params.pR, pl_all, az_R_hat, el_R_hat, az_T_hat, el_T_hat, ...
                    tau_hat, pR_hat, pl_hat);
                nmse_all(mc, 1, :) = nmse_vals;
                norm_sq_all(mc, :) = norms_sq;

                % MUSIC-LSPS
                if isfield(params, 'run_music') && params.run_music
                    [az_R_mu,el_R_mu,az_T_mu,el_T_mu,tau_mu,pR_mu,pl_mu] = ...
                        music_lsps(AR_hat, BT_hat, C_hat, params, pl_all, ...
                                   tau_true, az_R, el_R, az_T, el_T);
                    [nmse_mu, ~] = compute_nmse(az_R, el_R, az_T, el_T, tau_true, ...
                        params.pR, pl_all, az_R_mu, el_R_mu, az_T_mu, el_T_mu, ...
                        tau_mu, pR_mu, pl_mu);
                    nmse_all(mc, 2, :) = nmse_mu;
                end

                % PUDD
                if isfield(params, 'run_pudd') && params.run_pudd
                    [az_R_pu,el_R_pu,az_T_pu,el_T_pu,tau_pu,pR_pu,pl_pu] = ...
                        pudd_baseline(AR_hat, BT_hat, C_hat, params, pl_all, ...
                                      tau_true, az_R, el_R, az_T, el_T);
                    [nmse_pu, ~] = compute_nmse(az_R, el_R, az_T, el_T, tau_true, ...
                        params.pR, pl_all, az_R_pu, el_R_pu, az_T_pu, el_T_pu, ...
                        tau_pu, pR_pu, pl_pu);
                    nmse_all(mc, 3, :) = nmse_pu;
                end

                % CRB
                if isfield(params, 'run_crb') && params.run_crb
                    [CRB_p, CRB_pR_mat, CRB_pl_arr] = compute_crb(params, ...
                        AR_true, BT_true, C_true, pl_all, params.SNR_dB);
                    [crb_vec, crb_norms_sq] = extract_crb_nmse_local(CRB_p, CRB_pR_mat, CRB_pl_arr, ...
                        az_R, el_R, az_T, el_T, tau_true, params.pR, pl_all);
                    crb_all(mc, :) = crb_vec;
                    crb_norm_sq_all(mc, :) = crb_norms_sq;
                end
                n_valid = n_valid + 1;
            catch
            end
        end

        if n_valid == 0, n_valid = 1; end
        fprintf('    %d/%d valid trials\n', n_valid, MC);

        for mi_ = 1:length(methods)
            for fi_ = 1:length(fields)
                % Median MSE
                med_mse = nanmedian(nmse_all(:, mi_, fi_));
                % Convert to Median RMSE
                results.(methods{mi_}).(fields{fi_})(ni) = sqrt(med_mse);
            end
        end
        for fi_ = 1:length(fields)
            med_crb = nanmedian(crb_all(:, fi_));
            results.CRB.(fields{fi_})(ni) = sqrt(med_crb);
        end
    end
    results.N_vec = N_vec; results.params = params;
end

%% Local helpers

function [crb_vec, norms_sq] = extract_crb_nmse_local(CRB_p,CRB_pR,CRB_pl,az_R,el_R,az_T,el_T,tau,pR,pl)
    L=length(az_R);
    dim = size(CRB_p,1);
    if dim < 5*L
        crb_vec = zeros(1,7); norms_sq = zeros(1,7); return;
    end
    crb_vec=[sum(diag(CRB_p(1:L,1:L))), ...
             sum(diag(CRB_p(L+1:2*L,L+1:2*L))), ...
             sum(diag(CRB_p(2*L+1:3*L,2*L+1:3*L))), ...
             sum(diag(CRB_p(3*L+1:4*L,3*L+1:4*L))), ...
             sum(diag(CRB_p(4*L+1:5*L,4*L+1:5*L))), ...
             trace(CRB_pR), ...
             sum(arrayfun(@(l)trace(CRB_pl(:,:,l)),1:L))];
    norms_sq = [norm(az_R)^2, norm(el_R)^2, norm(az_T)^2, norm(el_T)^2, norm(tau)^2, norm(pR)^2, sum(sum(pl.^2))];
end

function perm = match_paths_combined(az_est, el_est, tau_est, az_true, el_true, tau_true)
    L = length(az_true);
    if L == 1, perm = 1; return; end
    cost_az  = abs(bsxfun(@minus, az_est(:), az_true(:)'));
    cost_el  = abs(bsxfun(@minus, el_est(:), el_true(:)'));
    cost_tau = abs(bsxfun(@minus, tau_est(:), tau_true(:)'));
    range_az  = max(eps, max(cost_az(:)));
    range_el  = max(eps, max(cost_el(:)));
    range_tau = max(eps, max(cost_tau(:)));
    cost = cost_az/range_az + cost_el/range_el + cost_tau/range_tau;
    perm = zeros(1, L); used = false(1, L);
    for i = 1:L
        c = cost(i,:); c(used) = inf;
        [~, j] = min(c); perm(i) = j; used(j) = true;
    end
end
