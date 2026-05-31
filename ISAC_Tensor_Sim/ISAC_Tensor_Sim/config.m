function params = config()
% =========================================================================
% config.m
% =========================================================================
% * Mathematical Background:
%   This function initializes all system and simulation parameters used in
%   the tensor-based ISAC simulation. It defines the mmWave MIMO-OFDM
%   system geometry, carrier frequency, bandwidth, and tensor dimensions
%   (F, T, K) as described in the paper.
%
% * Inputs:
%   None
%
% * Outputs:
%   params - A structure containing all simulation parameters:
%            - fc, c, fs, lambda, d (Physical parameters)
%            - NTy, NTz, NRy, NRz, NT, NR (Antenna arrays)
%            - MT, MR (RF chains)
%            - pT, pR (BS and UT positions)
%            - max_iter, eps_tol (ALS parameters)
%            - MC (Monte Carlo trials)
%
% * MATLAB Implementation:
%   A simple struct assignment.
%
% * Complexity Analysis:
%   O(1) time complexity and O(1) space complexity.
% =========================================================================

    % Physical parameters
    params.fc    = 30e9;          % Carrier frequency [Hz]
    params.c     = 3e8;           % Speed of light [m/s]
    params.fs    = 0.32e9;        % Sampling frequency [Hz]
    params.lambda = params.c / params.fc;  % Wavelength [m]
    params.d = params.lambda / 2;

    % Array sizes
    params.NTy = 7;   params.NTz = 7;
    params.NRy = 7;   params.NRz = 7;
    params.NT  = params.NTy * params.NTz;
    params.NR  = params.NRy * params.NRz;

    % RF chains
    params.MT  = 7;
    params.MR  = 1;

    % BS position (known anchor)
    params.pT  = [0; 0; 10];   % [m], 3x1

    % UT position (unknown, to be estimated)
    params.pR  = [10; 10; 0];  % [m], 3x1

    % SP search space size around the midpoint
    params.SP_range = 10; % meters

    % OFDM subcarriers total
    params.K_bar = 128;

    % Monte Carlo and algorithm settings
    params.MC  = 100; % Default MC trials
    
    % Toggle baseline algorithms
    params.run_music = false; % Set to false to save time
    params.run_pudd  = false; % Set to false to save time
    params.run_crb   = false; % Set to false to save time
    params.rng_seed = 42;
    params.Ns  = 1024;   % Number of 1D search points for ToA
    params.max_iter = 2000;
    params.eps_tol  = 1e-6;
end
