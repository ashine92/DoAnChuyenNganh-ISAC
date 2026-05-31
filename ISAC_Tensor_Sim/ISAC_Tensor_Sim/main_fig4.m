% =========================================================================
% main_fig4.m
% =========================================================================
% * Mathematical Background:
%   This script reproduces Figure 4 of the paper, evaluating the NMSE
%   (Normalized Mean Square Error) of estimated channel parameters (AoA,
%   AoD) and target localization (UT, SPs) versus Signal-to-Noise Ratio 
%   (SNR) ranging from 0 to 30 dB.
%
% * Inputs:
%   Parameters from config.m
%
% * Outputs:
%   Figure4_NMSE_vs_SNR.png - Plot matching the paper's Figure 4.
%
% * MATLAB Implementation:
%   Loops over SNR_dB_vec = 0:5:30. At each SNR, runs Monte Carlo trials.
%
% * Complexity Analysis:
%   O(N_SNR * MC * Complexity(run_monte_carlo))
% =========================================================================

clc; clear; close all;
addpath(genpath(pwd));

params = config();

% Override parameters specific to Figure 4
params.MC = 500; % As in paper for smooth monotonic curves
params.K = 10;   % Reverted to K=10, since ToA phase wrapping is now fixed
params.F = 49;   % F=NR for orthogonal DFT
params.T = 49;   % T=NT for orthogonal DFT
params.L = 3;
params.SNR_dB_vec = 0:5:30;

fprintf('--- Running Figure 4 (NMSE vs SNR) ---\n');
results4 = run_monte_carlo(params);

plot_fig4(results4, params);
fprintf('Figure 4 done.\n\n');
