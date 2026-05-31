% =========================================================================
% main_fig5.m
% =========================================================================
% * Mathematical Background:
%   This script reproduces Figure 5 of the paper, evaluating the NMSE
%   of estimated channel parameters and target localization versus the
%   number of selected subcarriers K (5 to 11).
%
% * Inputs:
%   Parameters from config.m
%
% * Outputs:
%   Figure5_NMSE_vs_K.png - Plot matching the paper's Figure 5.
%
% * MATLAB Implementation:
%   Loops over K_vec = 5:11 at a fixed SNR (20 dB), running MC trials.
%
% * Complexity Analysis:
%   O(N_K * MC * Complexity(run_monte_carlo))
% =========================================================================

clc; clear; close all;
addpath(genpath(pwd));

params = config();

% Override parameters specific to Figure 5
params.MC = 500; % As in paper for smooth monotonic curves
params.F = 49;  
params.T = 49;  
params.L = 3;  
params.SNR_dB = 20;
params.K_vec = 10:10:50;

fprintf('--- Running Figure 5 (NMSE vs K) ---\n');
results5 = run_monte_carlo_vs_K(params);

plot_fig5(results5, params);
fprintf('Figure 5 done.\n\n');
