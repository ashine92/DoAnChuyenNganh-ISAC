% =========================================================================
% quick_test.m
% =========================================================================
% Description:
%   Quick single-trial test to verify the full pipeline works end-to-end
%   before running the full Monte Carlo simulation.
%   Runtime: ~30-120 seconds.
%
% Run this FIRST before main.m to catch any issues.
% =========================================================================

clc; clear; close all;
addpath(genpath(pwd));

fprintf('=== Quick Pipeline Test ===\n');
rng(42);

%% System Parameters
params.fc    = 30e9;
params.c     = 3e8;
params.fs    = 0.32e9;
params.K_bar = 128;
params.lambda = params.c / params.fc;
params.d      = params.lambda / 4;
params.NTy = 7; params.NTz = 7;
params.NRy = 7; params.NRz = 7;
params.NT  = params.NTy * params.NTz;
params.NR  = params.NRy * params.NRz;
params.MT  = 7;  params.MR = 1;
params.pT  = [0; 0; 4*params.lambda];
params.pR  = [4*params.lambda; 4*params.lambda; 0];
params.SP_range = 4 * params.lambda;
params.F = 10; params.T = 10;   % Small for quick test
params.K = 5;  params.L = 2;
params.SNR_dB   = 20;
params.max_iter = 100;
params.eps_tol  = 1e-8;
params.Ns       = 256;
params.rng_seed = 42;

fprintf('Parameters: NT=%d, NR=%d, F=%d, T=%d, K=%d, L=%d, SNR=%d dB\n',...
    params.NT, params.NR, params.F, params.T, params.K, params.L, params.SNR_dB);

%% Step 1: Generate channel
fprintf('\n[1] Generating near-field channel... ');
tic;
[Hk, alpha, tau_true, az_R, el_R, az_T, el_T, pl_all, d_cR, d_cT] = generate_channel(params);
fprintf('Done (%.2f s)\n', toc);
fprintf('    tau_true = [%s] ns\n', num2str(tau_true'*1e9, '%.2f '));
fprintf('    az_R = [%s] deg\n', num2str(az_R'*180/pi, '%.1f '));
fprintf('    el_R = [%s] deg\n', num2str(el_R'*180/pi, '%.1f '));

%% Step 2: Verify near-field array response
fprintf('\n[2] Testing near_field_array_response... ');
[aR_test, bT_test] = near_field_array_response(params, pl_all(:,1));
fprintf('Done. |aR|=%.4f (expect ~%.4f)\n', norm(aR_test), sqrt(params.NR));

%% Step 3: Construct received tensor
fprintf('\n[3] Constructing received tensor... ');
tic;
[Y, W, F_mat, c_true, AR, BT] = construct_tensor(params, Hk, alpha, tau_true);
fprintf('Done (%.2f s)\n', toc);
fprintf('    Y size: %dx%dx%d, ||Y||_F = %.4f\n', size(Y,1), size(Y,2), size(Y,3), norm(Y(:)));

%% Step 4: MDL rank estimation
fprintf('\n[4] MDL rank estimation... ');
L_est = estimate_rank_mdl(Y, 5);
fprintf('Done. L_est = %d (true L = %d)\n', L_est, params.L);

%% Step 5: CP-ALS decomposition
fprintf('\n[5] CP-ALS decomposition (L=%d)... ', params.L);
tic;
[A_hat, B_hat, C_hat, iter_count, fit_history] = cp_als(Y, params.L, params);
fprintf('Done (%.2f s, %d iters, final fit = %.6f)\n', toc, iter_count, fit_history(end));

% Recover AR_hat, BT_hat
AR_hat = pinv(W') * A_hat;
BT_hat = pinv(F_mat') * B_hat;
fprintf('    ||AR_hat - AR_true||_F / ||AR_true||_F = %.4f\n', ...
    norm(AR_hat - AR_true,'fro') / norm(AR_true,'fro'));

%% Step 6: ToA estimation
fprintf('\n[6] ToA estimation... ');
tic;
tau_hat = estimate_toa(C_hat, params);
fprintf('Done (%.2f s)\n', toc);
fprintf('    tau_true = [%s] ns\n', num2str(sort(tau_true)'*1e9, '%.2f '));
fprintf('    tau_hat  = [%s] ns\n', num2str(sort(tau_hat)'*1e9, '%.2f '));
nmse_tau = norm(sort(tau_true)-sort(tau_hat))^2 / norm(tau_true)^2;
fprintf('    NMSE(tau) = %.4e\n', nmse_tau);

%% Step 7: Angle estimation (AoA)
fprintf('\n[7] AoA estimation... ');
tic;
[az_R_hat, el_R_hat] = estimate_angles(AR_hat, params);
fprintf('Done (%.2f s)\n', toc);
[az_R_s, ia] = sort(az_R); el_R_s = el_R(ia);
[az_R_hat_s, ib] = sort(az_R_hat); el_R_hat_s = el_R_hat(ib);
fprintf('    az_R true = [%s] deg\n', num2str(az_R_s'*180/pi, '%.2f '));
fprintf('    az_R hat  = [%s] deg\n', num2str(az_R_hat_s'*180/pi, '%.2f '));
fprintf('    el_R true = [%s] deg\n', num2str(el_R_s'*180/pi, '%.2f '));
fprintf('    el_R hat  = [%s] deg\n', num2str(el_R_hat_s'*180/pi, '%.2f '));

%% Step 8: Angle estimation (AoD)
fprintf('\n[8] AoD estimation... ');
params_T = params; params_T.NRy=params.NTy; params_T.NRz=params.NTz; params_T.NR=params.NT;
[az_T_hat, el_T_hat] = estimate_angles(BT_hat, params_T);
fprintf('Done\n');
[az_T_s, ic] = sort(az_T); el_T_s = el_T(ic);
[az_T_hat_s, id] = sort(az_T_hat); el_T_hat_s = el_T_hat(id);
fprintf('    az_T true = [%s] deg\n', num2str(az_T_s'*180/pi, '%.2f '));
fprintf('    az_T hat  = [%s] deg\n', num2str(az_T_hat_s'*180/pi, '%.2f '));

%% Step 9: UT Localization
fprintf('\n[9] UT localization... ');
pR_hat = ut_localization(sort(tau_hat), az_R_hat_s, el_R_hat_s, az_T_hat_s, el_T_hat_s, params);
fprintf('Done\n');
fprintf('    pR true = [%.4f, %.4f, %.4f] m\n', params.pR');
fprintf('    pR hat  = [%.4f, %.4f, %.4f] m\n', pR_hat');
nmse_pR = norm(params.pR - pR_hat)^2 / norm(params.pR)^2;
fprintf('    NMSE(pR) = %.4e\n', nmse_pR);

%% Step 10: SP Localization
fprintf('\n[10] SP localization... ');
pl_hat = sp_localization(sort(tau_hat), az_R_hat_s, el_R_hat_s, az_T_hat_s, el_T_hat_s, pR_hat, params);
fprintf('Done\n');
for l = 1:params.L
    fprintf('    SP%d true = [%.4f, %.4f, %.4f] m\n', l, pl_all(:,l)');
    fprintf('    SP%d hat  = [%.4f, %.4f, %.4f] m\n', l, pl_hat(:,l)');
end
nmse_pl = sum(sum((pl_all - pl_hat).^2, 1)) / sum(sum(pl_all.^2, 1));
fprintf('    NMSE(pl) = %.4e\n', nmse_pl);

%% Step 11: CRB
fprintf('\n[11] Computing CRB... ');
tic;
AR_true_all = zeros(params.NR, params.L);
BT_true_all = zeros(params.NT, params.L);
for l = 1:params.L
    [AR_true_all(:,l), BT_true_all(:,l)] = near_field_array_response(params, pl_all(:,l));
end
k_idx = round(linspace(1, params.K_bar, params.K))';
C_true_all = zeros(params.K, params.L);
for l = 1:params.L
    C_true_all(:,l) = alpha(l)*exp(-1j*2*pi*tau_true(l).*params.fs.*k_idx/params.K_bar);
end
[CRB_params, CRB_pR, CRB_pl_arr] = compute_crb(params, AR_true_all, BT_true_all, C_true_all, pl_all, params.SNR_dB);
fprintf('Done (%.2f s)\n', toc);
fprintf('    CRB trace for [az_R, el_R, az_T, el_T] = [%.2e, %.2e, %.2e, %.2e]\n', ...
    trace(CRB_params(1:params.L,1:params.L)), ...
    trace(CRB_params(params.L+1:2*params.L,params.L+1:2*params.L)), ...
    trace(CRB_params(2*params.L+1:3*params.L,2*params.L+1:3*params.L)), ...
    trace(CRB_params(3*params.L+1:4*params.L,3*params.L+1:4*params.L)));
fprintf('    CRB trace for pR = %.2e\n', trace(CRB_pR));

%% Step 12: 3D visualization
fprintf('\n[12] 3D visualization... ');
figure('Name','Quick Test - 3D Localization');
hold on; grid on; box on; view(3);
plot3(params.pT(1), params.pT(2), params.pT(3), 'rs','MarkerSize',12,'MarkerFaceColor','r','DisplayName','BS');
plot3(params.pR(1), params.pR(2), params.pR(3), 'ko','MarkerSize',10,'MarkerFaceColor','k','DisplayName','UT true');
plot3(pR_hat(1), pR_hat(2), pR_hat(3), 'g^','MarkerSize',10,'MarkerFaceColor','g','DisplayName','UT est');
plot3([params.pR(1),pR_hat(1)],[params.pR(2),pR_hat(2)],[params.pR(3),pR_hat(3)],'g--','LineWidth',1,'HandleVisibility','off');
for l = 1:params.L
    plot3(pl_all(1,l),pl_all(2,l),pl_all(3,l),'bs','MarkerSize',8,'MarkerFaceColor','b','DisplayName',sprintf('SP%d true',l));
    plot3(pl_hat(1,l),pl_hat(2,l),pl_hat(3,l),'ms','MarkerSize',8,'MarkerFaceColor','m','DisplayName',sprintf('SP%d est',l));
    plot3([pl_all(1,l),pl_hat(1,l)],[pl_all(2,l),pl_hat(2,l)],[pl_all(3,l),pl_hat(3,l)],'m--','LineWidth',1,'HandleVisibility','off');
end
legend('Location','best');
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(sprintf('Quick Test: SNR=%d dB, L=%d', params.SNR_dB, params.L));
fprintf('Done\n');

%% Summary
fprintf('\n=== Pipeline Test Summary ===\n');
fprintf('  NMSE(tau) = %.4e\n', nmse_tau);
fprintf('  NMSE(az_R) = %.4e\n', norm(sort(az_R)-az_R_hat_s)^2/norm(az_R)^2);
fprintf('  NMSE(pR)   = %.4e\n', nmse_pR);
fprintf('  NMSE(pl)   = %.4e\n', nmse_pl);
fprintf('\nAll steps passed! Ready to run main.m\n');
fprintf('(Reduce params.MC in main.m for faster testing, e.g. MC=20)\n');
