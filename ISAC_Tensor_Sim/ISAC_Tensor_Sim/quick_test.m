% =========================================================================
% quick_test.m
% =========================================================================
% Quick single-trial pipeline test. Run this BEFORE main.m.
% Expected runtime: ~30-90 seconds.
%
% Checks every stage: channel gen → tensor → CP-ALS → ToA → angles →
% localization → CRB → 3D plot.
% =========================================================================

clc; clear; close all;
addpath(genpath(pwd));

fprintf('=== Quick Pipeline Test ===\n');
rng(42, 'twister');

%% -----------------------------------------------------------------------
% SYSTEM PARAMETERS (small F,T,K for speed)
% -----------------------------------------------------------------------
params.fc     = 30e9;
params.c      = 3e8;
params.fs     = 0.32e9;
params.K_bar  = 128;
params.lambda = params.c / params.fc;
params.d      = params.lambda / 4;
params.NTy = 7; params.NTz = 7;
params.NRy = 7; params.NRz = 7;
params.NT  = params.NTy * params.NTz;   % 49
params.NR  = params.NRy * params.NRz;   % 49
params.MT  = 7;  params.MR = 1;
params.pT  = [0; 0; 4*params.lambda];
params.pR  = [4*params.lambda; 4*params.lambda; 0];
params.SP_range   = 4 * params.lambda;
params.F          = 10;   % small for quick test
params.T          = 10;
params.K          = 6;
params.L          = 2;
params.SNR_dB     = 20;
params.max_iter   = 200;
params.eps_tol    = 1e-8;
params.Ns         = 512;
params.rng_seed   = 42;

fprintf('NT=%d, NR=%d  |  F=%d, T=%d, K=%d, L=%d, SNR=%d dB\n', ...
    params.NT, params.NR, params.F, params.T, params.K, params.L, params.SNR_dB);

%% -----------------------------------------------------------------------
% STEP 1 – Channel generation
% -----------------------------------------------------------------------
fprintf('\n[1] generate_channel ... ');
tic;
[Hk, alpha, tau_true, az_R, el_R, az_T, el_T, pl_all, d_cR, d_cT] = ...
    generate_channel(params);
fprintf('%.2f s\n', toc);
fprintf('    tau_true = [%s] ns\n', num2str(tau_true'*1e9, '%.3f '));
fprintf('    az_R_true = [%s] deg\n', num2str(az_R'*180/pi, '%.2f '));
fprintf('    el_R_true = [%s] deg\n', num2str(el_R'*180/pi, '%.2f '));

%% -----------------------------------------------------------------------
% STEP 2 – Array responses (build AR_true / BT_true FIRST)
% -----------------------------------------------------------------------
fprintf('\n[2] near_field_array_response ... ');
AR_true = zeros(params.NR, params.L);
BT_true = zeros(params.NT, params.L);
for l = 1:params.L
    [AR_true(:,l), BT_true(:,l)] = near_field_array_response(params, pl_all(:,l));
end
fprintf('Done.  |aR_1| = %.4f (expect sqrt(NR)=%.4f)\n', ...
    norm(AR_true(:,1)), sqrt(params.NR));

%% -----------------------------------------------------------------------
% STEP 3 – Construct received tensor (pass pl_all so AR,BT are returned)
% -----------------------------------------------------------------------
fprintf('\n[3] construct_tensor ... ');
tic;
[Y, W, F_mat, c_true, AR_ret, BT_ret] = construct_tensor(params, Hk, alpha, tau_true, pl_all);
fprintf('%.2f s  |  Y size: %dx%dx%d,  ||Y||_F = %.4f\n', ...
    toc, size(Y,1), size(Y,2), size(Y,3), norm(Y(:)));

% Verify AR from construct_tensor matches the one we built above
if ~isempty(AR_ret)
    ar_err = norm(AR_ret - AR_true,'fro') / norm(AR_true,'fro');
    fprintf('    AR consistency check: %.2e (expect ~0)\n', ar_err);
end

% True C factor for later CRB
k_indices = round(linspace(1, params.K_bar, params.K))';
C_true = zeros(params.K, params.L);
for l = 1:params.L
    C_true(:,l) = alpha(l)*exp(-1j*2*pi*tau_true(l).*params.fs.*k_indices/params.K_bar);
end

%% -----------------------------------------------------------------------
% STEP 4 – MDL rank estimation
% -----------------------------------------------------------------------
fprintf('\n[4] MDL rank estimation ... ');
L_est = estimate_rank_mdl(Y, 5);
fprintf('L_est = %d  (true L = %d)\n', L_est, params.L);

%% -----------------------------------------------------------------------
% STEP 5 – CP-ALS decomposition
% -----------------------------------------------------------------------
fprintf('\n[5] CP-ALS (L=%d) ... ', params.L);
tic;
[A_hat, B_hat, C_hat, iter_count, fit_history] = cp_als(Y, params.L, params);
fprintf('%.2f s  |  %d iters, final fit = %.6f\n', toc, iter_count, fit_history(end));
if fit_history(end) > 1
    fprintf('    WARNING: fit > 1 — ALS did not converge well at this SNR/size.\n');
end

% Recover AR_hat, BT_hat via pseudo-inverse of DFT slices (Eq. 31)
% Recover AR_hat, BT_hat via Eq.(31): A = W'*AR  →  AR ≈ W*A  (pinv(W')=W for partial DFT)
AR_hat = W * A_hat;     % NR×L
BT_hat = F_mat * B_hat; % NT×L

fprintf('    ||AR_hat - AR_true||_F/||AR_true||_F = %.4f\n', ...
    norm(AR_hat - AR_true,'fro') / norm(AR_true,'fro'));

%% -----------------------------------------------------------------------
% STEP 6 – ToA estimation
% -----------------------------------------------------------------------
fprintf('\n[6] estimate_toa ... ');
tic;
tau_hat = estimate_toa(C_hat, params);
fprintf('%.2f s\n', toc);
[tau_s, ~] = sort(tau_true);
[tau_h, ~] = sort(tau_hat);
fprintf('    tau_true = [%s] ns\n', num2str(tau_s'*1e9, '%.3f '));
fprintf('    tau_hat  = [%s] ns\n', num2str(tau_h'*1e9, '%.3f '));
nmse_tau = norm(tau_s - tau_h)^2 / (norm(tau_s)^2 + eps);
fprintf('    NMSE(tau) = %.4e\n', nmse_tau);

% Path matching by ToA
[tau_hat, perm] = match_paths_qt(tau_hat, tau_true);
AR_hat = AR_hat(:, perm);
BT_hat = BT_hat(:, perm);

%% -----------------------------------------------------------------------
% STEP 7 – AoA estimation
% -----------------------------------------------------------------------
fprintf('\n[7] estimate_angles (AoA) ... ');
[az_R_hat, el_R_hat] = estimate_angles(AR_hat, params);
[az_R_sorted, ia] = sort(az_R);  el_R_sorted = el_R(ia);
[az_R_hat_s,  ib] = sort(real(az_R_hat));  el_R_hat_s = real(el_R_hat(ib));
fprintf('Done\n');
fprintf('    az_R true = [%s] deg\n', num2str(az_R_sorted'*180/pi,  '%.2f '));
fprintf('    az_R hat  = [%s] deg\n', num2str(az_R_hat_s'*180/pi,   '%.2f '));
fprintf('    el_R true = [%s] deg\n', num2str(el_R_sorted'*180/pi,  '%.2f '));
fprintf('    el_R hat  = [%s] deg\n', num2str(el_R_hat_s'*180/pi,   '%.2f '));

%% -----------------------------------------------------------------------
% STEP 8 – AoD estimation
% -----------------------------------------------------------------------
fprintf('\n[8] estimate_angles (AoD) ... ');
params_T = params;
params_T.NRy = params.NTy;  params_T.NRz = params.NTz;  params_T.NR = params.NT;
[az_T_hat, el_T_hat] = estimate_angles(BT_hat, params_T);
[az_T_sorted, ic] = sort(az_T);  el_T_sorted = el_T(ic);
[az_T_hat_s,  id] = sort(real(az_T_hat));  el_T_hat_s = real(el_T_hat(id));
fprintf('Done\n');
fprintf('    az_T true = [%s] deg\n', num2str(az_T_sorted'*180/pi,  '%.2f '));
fprintf('    az_T hat  = [%s] deg\n', num2str(az_T_hat_s'*180/pi,   '%.2f '));

%% -----------------------------------------------------------------------
% STEP 9 – UT Localization
% -----------------------------------------------------------------------
fprintf('\n[9] ut_localization ... ');
pR_hat = ut_localization(tau_hat, az_R_hat_s, el_R_hat_s, az_T_hat_s, el_T_hat_s, params);
nmse_pR = norm(params.pR - pR_hat)^2 / (norm(params.pR)^2 + eps);
fprintf('Done\n');
fprintf('    pR true = [%.4f, %.4f, %.4f] m\n', params.pR');
fprintf('    pR hat  = [%.4f, %.4f, %.4f] m\n', pR_hat');
fprintf('    NMSE(pR) = %.4e\n', nmse_pR);

%% -----------------------------------------------------------------------
% STEP 10 – SP Localization
% -----------------------------------------------------------------------
fprintf('\n[10] sp_localization ... ');
pl_hat = sp_localization(tau_hat, az_R_hat_s, el_R_hat_s, az_T_hat_s, el_T_hat_s, pR_hat, params);
nmse_pl = sum(sum((pl_all - pl_hat).^2,1)) / (sum(sum(pl_all.^2,1)) + eps);
fprintf('Done\n');
for l = 1:params.L
    fprintf('    SP%d true = [%.4f, %.4f, %.4f] m\n', l, pl_all(:,l)');
    fprintf('    SP%d hat  = [%.4f, %.4f, %.4f] m\n', l, pl_hat(:,l)');
end
fprintf('    NMSE(pl) = %.4e\n', nmse_pl);

%% -----------------------------------------------------------------------
% STEP 11 – CRB
% -----------------------------------------------------------------------
fprintf('\n[11] compute_crb ... ');
tic;
[CRB_params, CRB_pR_mat, CRB_pl_arr] = compute_crb(params, AR_true, BT_true, C_true, pl_all, params.SNR_dB);
fprintf('%.2f s\n', toc);
fprintf('    CRB(az_R) trace = %.2e\n', trace(CRB_params(1:params.L, 1:params.L)));
fprintf('    CRB(pR)   trace = %.2e\n', trace(CRB_pR_mat));

%% -----------------------------------------------------------------------
% STEP 12 – 3D visualization
% -----------------------------------------------------------------------
figure('Name','Quick Test – 3D Localization','Position',[100 100 600 500]);
hold on; grid on; box on; view([-35 30]);

plot3(params.pT(1), params.pT(2), params.pT(3), 'rs','MarkerSize',12,...
    'MarkerFaceColor','r','DisplayName','BS');
plot3(params.pR(1), params.pR(2), params.pR(3), 'ko','MarkerSize',10,...
    'MarkerFaceColor','k','DisplayName','UT true');
plot3(pR_hat(1), pR_hat(2), pR_hat(3), 'g^','MarkerSize',10,...
    'MarkerFaceColor','g','DisplayName','UT est');
plot3([params.pR(1),pR_hat(1)],[params.pR(2),pR_hat(2)],[params.pR(3),pR_hat(3)],...
    'g--','LineWidth',1.5,'HandleVisibility','off');

for l = 1:params.L
    plot3(pl_all(1,l),pl_all(2,l),pl_all(3,l),'bs','MarkerSize',9,...
        'MarkerFaceColor','b','DisplayName',sprintf('SP%d true',l));
    plot3(pl_hat(1,l),pl_hat(2,l),pl_hat(3,l),'ms','MarkerSize',9,...
        'MarkerFaceColor','m','DisplayName',sprintf('SP%d est',l));
    plot3([pl_all(1,l),pl_hat(1,l)],[pl_all(2,l),pl_hat(2,l)],...
        [pl_all(3,l),pl_hat(3,l)],'m--','LineWidth',1.5,'HandleVisibility','off');
end
legend('Location','best','FontSize',9);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(sprintf('Quick Test: SNR=%d dB, L=%d, F=T=%d, K=%d',...
    params.SNR_dB, params.L, params.F, params.K));

%% -----------------------------------------------------------------------
% SUMMARY
% -----------------------------------------------------------------------
fprintf('\n======= PIPELINE TEST SUMMARY =======\n');
fprintf('  NMSE(tau)  = %.4e\n', nmse_tau);
nmse_azR = norm(az_R_sorted - az_R_hat_s)^2 / (norm(az_R_sorted)^2+eps);
fprintf('  NMSE(az_R) = %.4e\n', nmse_azR);
fprintf('  NMSE(pR)   = %.4e\n', nmse_pR);
fprintf('  NMSE(pl)   = %.4e\n', nmse_pl);
fprintf('=====================================\n');
fprintf('  ALS fit = %.4f  (good if < 0.5 at SNR=20 dB)\n', fit_history(end));
fprintf('\n  All steps completed. Ready to run main.m\n');
fprintf('  (Set params.MC = 20 in main.m for a fast preview)\n');


%% ---- Local helper: greedy path matching --------------------------------
function [x_matched, perm] = match_paths_qt(x_est, x_true)
    L = length(x_true);
    if L == 1, perm = 1; x_matched = x_est; return; end
    cost = abs(bsxfun(@minus, real(x_est(:)), real(x_true(:))'));
    perm = zeros(1,L);  used = false(1,L);
    for i = 1:L
        [~,j] = min(cost(i,:) + 1e9*used);
        perm(i) = j;  used(j) = true;
    end
    x_matched = x_est(perm);
end
