% =========================================================================
% run_fig4_quick.m — Quick Figure 4 with MC=30 for initial validation
% =========================================================================
clc; clear; close all;
addpath(genpath(pwd));

params = config();
params.F = 49;  % F=NR to ensure W'W=I (orthogonal DFT)
params.T = 49;  % T=NT to ensure F'F=I
params.K = 10;  % Restored to 10 as in paper
params.L = 3;
params.MC = 300;
params.SNR_dB_vec = 0:5:30;
params.run_music = false;
params.run_pudd = false;
params.run_crb = true;

fprintf('--- Running Figure 4 (NMSE vs SNR, MC=%d) ---\n', params.MC);
tic;
results4 = run_monte_carlo(params);
t = toc;
fprintf('Total time: %.1f sec\n', t);

% Print results table
fprintf('\n=== NMSE vs SNR Results ===\n');
fprintf('%5s', 'SNR');
fields = {'az_R','el_R','az_T','el_T','tau','pR','pl'};
for f = fields, fprintf('%12s', f{1}); end
fprintf('\n');

for si = 1:length(results4.SNR_vec)
    fprintf('%5d', results4.SNR_vec(si));
    for f = fields
        fprintf('%12.4e', results4.Proposed.(f{1})(si));
    end
    fprintf('\n');
end

fprintf('\nCRB:\n');
for si = 1:length(results4.SNR_vec)
    fprintf('%5d', results4.SNR_vec(si));
    for f = fields
        fprintf('%12.4e', results4.CRB.(f{1})(si));
    end
    fprintf('\n');
end

% Self-check: NMSE should decrease with SNR
fprintf('\n=== Self-Check ===\n');
for f = fields
    vals = results4.Proposed.(f{1});
    is_decreasing = all(diff(vals(3:end)) <= 0.01);  % Allow small fluctuation
    if is_decreasing
        fprintf('CHECK %s decreases with SNR: PASS\n', f{1});
    else
        fprintf('CHECK %s decreases with SNR: FAIL (values: %s)\n', f{1}, num2str(vals, '%.3e '));
    end
end

plot_fig4(results4, params);
