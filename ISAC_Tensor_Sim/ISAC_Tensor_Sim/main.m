% =========================================================================
% main.m
% =========================================================================
% Description:
%   Main simulation script for Near-Field Channel Parameter Estimation and
%   Localization for mmWave Massive MIMO-OFDM ISAC Systems via Tensor Analysis.
%   Reproduces Figures 4-8 from:
%   Jiang et al., "Near-Field Channel Parameter Estimation and Localization
%   for mmWave Massive MIMO-OFDM ISAC Systems via Tensor Analysis,"
%   Sensors 2025, 25, 5050.
%
% Inputs:  None (parameters set in params struct)
% Outputs: Figures 4-8 saved as PNG files
%
% Usage:   Run this script directly from MATLAB.
%          Ensure all .m files are on the MATLAB path.
%
% Reproducibility: rng(42) seed is set for Monte Carlo trials.
% Runtime: ~20-60 min for full MC=600; reduce MC for faster testing.
% =========================================================================

clc; clear; close all;
addpath(genpath(pwd));

%% -----------------------------------------------------------------------
% SYSTEM PARAMETERS (Table 1 / Section 7 of paper)
% -----------------------------------------------------------------------
params.fc    = 30e9;          % Carrier frequency [Hz]
params.c     = 3e8;           % Speed of light [m/s]
params.fs    = 0.32e9;        % Sampling frequency [Hz]
params.K_bar = 128;           % Total OFDM subcarriers
params.lambda = params.c / params.fc;  % Wavelength [m]
params.d     = params.lambda / 4;      % Antenna spacing [m]

% Array sizes
params.NTy = 7;   params.NTz = 7;
params.NRy = 7;   params.NRz = 7;
params.NT  = params.NTy * params.NTz;   % = 49
params.NR  = params.NRy * params.NRz;   % = 49

% RF chains
params.MT  = 7;   % < NT; chosen as NTy for DFT construction
params.MR  = 1;

% BS position (known anchor)
params.pT  = [0; 0; 4*params.lambda];   % [m], 3×1

% UT position (unknown, to be estimated)
params.pR  = [4*params.lambda; 4*params.lambda; 0];  % [m], 3×1

% SP search space: 4λ × 4λ × 4λ cube (relative coords used in generation)
params.SP_range = 4 * params.lambda;

% Noise
params.noise_type = 'AWGN';

% Monte Carlo
params.MC  = 600;    % Number of trials (reduce to 50 for quick test)
params.rng_seed = 42;

% ToA search
params.Ns  = 1024;   % Number of 1D search points for ToA

% ALS parameters
params.max_iter = 500;
params.eps_tol  = 1e-10;

% -----------------------------------------------------------------------
% EXPERIMENT SELECTION
% -----------------------------------------------------------------------
run_fig4 = true;   % NMSE vs SNR  (F=T=50, L=3, K=10)
run_fig5 = true;   % NMSE vs K    (F=T=50, L=3, SNR=20)
run_fig6 = true;   % NMSE vs SNR  (L=2,3 ; F=T=50,80)
run_fig7 = true;   % ToA/Loc NMSE vs SNR
run_fig8 = true;   % 3D visualization

fprintf('=== Near-Field ISAC Tensor Simulation ===\n');
fprintf('NT=%d, NR=%d, fc=%.0f GHz, fs=%.2f GHz\n', ...
    params.NT, params.NR, params.fc/1e9, params.fs/1e9);
fprintf('MC trials = %d\n\n', params.MC);

%% -----------------------------------------------------------------------
% FIGURE 4 : NMSE vs SNR  (F=T=50, K=10, L=3)
% -----------------------------------------------------------------------
if run_fig4
    fprintf('--- Running Figure 4 (NMSE vs SNR) ---\n');
    p4 = params;
    p4.F = 50;  p4.T = 50;  p4.K = 10;  p4.L = 3;
    p4.SNR_dB_vec = 0:5:30;
    results4 = run_monte_carlo(p4);
    plot_fig4(results4, p4);
    fprintf('Figure 4 done.\n\n');
end

%% -----------------------------------------------------------------------
% FIGURE 5 : NMSE vs K  (F=T=50, L=3, SNR=20 dB)
% -----------------------------------------------------------------------
if run_fig5
    fprintf('--- Running Figure 5 (NMSE vs K) ---\n');
    p5 = params;
    p5.F = 50;  p5.T = 50;  p5.L = 3;  p5.SNR_dB = 20;
    p5.K_vec = 5:11;
    results5 = run_monte_carlo_vs_K(p5);
    plot_fig5(results5, p5);
    fprintf('Figure 5 done.\n\n');
end

%% -----------------------------------------------------------------------
% FIGURES 6 & 7 : NMSE vs SNR for 4 cases (K=8)
% -----------------------------------------------------------------------
if run_fig6 || run_fig7
    fprintf('--- Running Figures 6&7 (4 cases) ---\n');
    p67 = params;
    p67.K = 8;
    p67.SNR_dB_vec = 0:5:30;
    cases = struct('L', {2,2,3,3}, 'F', {80,50,80,50}, 'T', {80,50,80,50}, ...
                   'label', {'L=2,F=T=80','L=2,F=T=50','L=3,F=T=80','L=3,F=T=50'});
    results67 = cell(1,4);
    for ci = 1:4
        fprintf('  Case %d: %s\n', ci, cases(ci).label);
        p67.L = cases(ci).L;
        p67.F = cases(ci).F;
        p67.T = cases(ci).T;
        results67{ci} = run_monte_carlo(p67);
    end
    if run_fig6
        plot_fig6(results67, cases, p67);
        fprintf('Figure 6 done.\n\n');
    end
    if run_fig7
        plot_fig7(results67, cases, p67);
        fprintf('Figure 7 done.\n\n');
    end
end

%% -----------------------------------------------------------------------
% FIGURE 8 : 3D Localization Visualization
% -----------------------------------------------------------------------
if run_fig8
    fprintf('--- Running Figure 8 (3D Localization) ---\n');
    p8 = params;
    p8.K = 8;  p8.SNR_dB = 5;
    cases8 = struct('L',{3,3,2,2},'F',{50,80,50,80},'T',{50,80,50,80});
    results8 = cell(1,4);
    for ci = 1:4
        p8.L = cases8(ci).L;  p8.F = cases8(ci).F;  p8.T = cases8(ci).T;
        results8{ci} = run_single_trial(p8);
    end
    plot_fig8(results8, cases8, p8);
    fprintf('Figure 8 done.\n\n');
end

fprintf('=== Simulation Complete ===\n');
