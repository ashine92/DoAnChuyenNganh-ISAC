% =========================================================================
% run_fig6_quick.m — Quick Figure 6 with MC=300 for validation
% =========================================================================
clc; clear; close all;
addpath(genpath(pwd));

params = config();
params.L = 3;
params.MC = 300;
params.SNR_dB = 10;
params.K = 10;
params.N_vec = 7:2:15; % N_y = N_z = N
params.run_music = false;
params.run_pudd = false;
params.run_crb = true;
fprintf('--- Running Figure 6 (RMSE vs N, MC=%d) ---\n', params.MC);
tic;
results6 = run_monte_carlo_vs_N(params);
t = toc;
fprintf('Total time: %.1f sec\n', t);

% Print results table
fprintf('\n=== RMSE vs N Results ===\n');
fprintf('%5s', 'N');
fields = {'az_R','el_R','az_T','el_T','pR','pl'};
for f = fields, fprintf('%12s', f{1}); end
fprintf('\n');

for ni = 1:length(results6.N_vec)
    N_val = results6.N_vec(ni)^2;
    fprintf('%5d', N_val);
    for f = fields
        fprintf('%12.4e', results6.Proposed.(f{1})(ni));
    end
    fprintf('\n');
end

% Self-check: NMSE should generally decrease with K
fprintf('\n=== Self-Check ===\n');
for f = fields
    vals = results6.Proposed.(f{1});
    trend = polyfit(results6.N_vec, log10(max(vals, 1e-12)), 1);
    if trend(1) < 0
        fprintf('CHECK %s decreasing trend with N: PASS (slope=%.4f)\n', f{1}, trend(1));
    else
        fprintf('CHECK %s decreasing trend with N: FAIL (slope=%.4f)\n', f{1}, trend(1));
    end
end

% Plot and save
plot_fig6(results6, params);
