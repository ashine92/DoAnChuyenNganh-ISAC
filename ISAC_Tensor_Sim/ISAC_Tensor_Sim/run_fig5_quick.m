% =========================================================================
% run_fig5_quick.m — Quick Figure 5 with MC=30 for initial validation
% =========================================================================
clc; clear; close all;
addpath(genpath(pwd));

params = config();
params.F = 49;
params.T = 49;
params.L = 3;
params.MC = 300;
params.SNR_dB = 15;
params.K_vec = 10:10:50;
params.run_music = false;
params.run_pudd = false;
params.run_crb = true;
fprintf('--- Running Figure 5 (NMSE vs K, MC=%d) ---\n', params.MC);
tic;
results5 = run_monte_carlo_vs_K(params);
t = toc;
fprintf('Total time: %.1f sec\n', t);

% Print results table
fprintf('\n=== NMSE vs K Results ===\n');
fprintf('%5s', 'K');
fields = {'az_R','el_R','az_T','el_T','pR','pl'};
for f = fields, fprintf('%12s', f{1}); end
fprintf('\n');

for ki = 1:length(results5.K_vec)
    fprintf('%5d', results5.K_vec(ki));
    for f = fields
        fprintf('%12.4e', results5.Proposed.(f{1})(ki));
    end
    fprintf('\n');
end

% Self-check: NMSE should generally decrease with K
fprintf('\n=== Self-Check ===\n');
for f = fields
    vals = results5.Proposed.(f{1});
    trend = polyfit(results5.K_vec, log10(max(vals, 1e-12)), 1);
    if trend(1) < 0
        fprintf('CHECK %s decreasing trend with K: PASS (slope=%.4f)\n', f{1}, trend(1));
    else
        fprintf('CHECK %s decreasing trend with K: FAIL (slope=%.4f)\n', f{1}, trend(1));
    end
end

plot_fig5(results5, params);
