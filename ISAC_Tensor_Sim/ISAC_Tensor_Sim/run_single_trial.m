function result = run_single_trial(params)
% =========================================================================
% run_single_trial.m
% =========================================================================
% Description:
%   Runs a single Monte Carlo trial at fixed SNR for 3D localization
%   visualization (Figure 8). Returns true and estimated positions.
%
% Inputs:
%   params - system parameter struct with SNR_dB set
%
% Outputs:
%   result - struct with true/estimated positions and params
% =========================================================================

    rng(params.rng_seed + 100, 'twister');

    L = params.L;

    % --- Generate channel ---
    [Hk, alpha, tau_true, az_R, el_R, az_T, el_T, pl_all, ~, ~] = ...
        generate_channel(params);

    % --- Construct tensor ---
    [Y, W, F_mat, ~, ~, ~] = construct_tensor(params, Hk, alpha, tau_true);

    % --- CP-ALS ---
    [A_hat, B_hat, C_hat, ~, ~] = cp_als(Y, L, params);
    AR_hat = pinv(W') * A_hat;
    BT_hat = pinv(F_mat') * B_hat;

    % --- ToA ---
    tau_hat = estimate_toa(C_hat, params);
    [tau_hat, perm] = match_paths_st(tau_hat, tau_true);
    AR_hat = AR_hat(:,perm); BT_hat = BT_hat(:,perm);

    % --- Angles ---
    params_R = params;
    [az_R_hat, el_R_hat] = estimate_angles(AR_hat, params_R);
    params_T = params; params_T.NRy=params.NTy; params_T.NRz=params.NTz; params_T.NR=params.NT;
    [az_T_hat, el_T_hat] = estimate_angles(BT_hat, params_T);
    [az_R_hat, p2] = match_paths_st(az_R_hat, az_R); el_R_hat=el_R_hat(p2);
    [az_T_hat, p3] = match_paths_st(az_T_hat, az_T); el_T_hat=el_T_hat(p3);

    % --- Localization ---
    pR_hat = ut_localization(tau_hat, az_R_hat, el_R_hat, az_T_hat, el_T_hat, params);
    pl_hat = sp_localization(tau_hat, az_R_hat, el_R_hat, az_T_hat, el_T_hat, pR_hat, params);

    % --- Store results ---
    result.pT        = params.pT;
    result.pR_true   = params.pR;
    result.pR_hat    = pR_hat;
    result.pl_true   = pl_all;
    result.pl_hat    = pl_hat;
    result.tau_true  = tau_true;
    result.tau_hat   = tau_hat;
    result.L         = L;
    result.params    = params;
end

function [x_matched, perm] = match_paths_st(x_est, x_true)
    L = length(x_true);
    if L==1, perm=1; x_matched=x_est; return; end
    cost = abs(bsxfun(@minus, real(x_est(:)), real(x_true(:))'));
    perm = zeros(1,L); used = false(1,L);
    for i=1:L
        [~,j] = min(cost(i,:)+1e9*used);
        perm(i)=j; used(j)=true;
    end
    x_matched = x_est(perm);
end
