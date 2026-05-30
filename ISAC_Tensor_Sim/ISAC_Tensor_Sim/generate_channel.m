function [Hk_all, alpha, tau, theta_az_R, theta_el_R, theta_az_T, theta_el_T, pl_all, d_c_R, d_c_T] = ...
    generate_channel(params)
% =========================================================================
% generate_channel.m
% =========================================================================
% Description:
%   Generates the near-field mmWave MIMO-OFDM channel for all K selected
%   subcarriers. Uses the NLoS-only model with spherical wavefront.
%
%   Time-domain impulse response (Eq. 11):
%     H(tau) = sum_l alpha_l * aR_l * bT_l^T * delta(tau - tau_l)
%
%   Frequency-domain at subcarrier k (Eq. 14):
%     Hk = sum_l alpha_l * exp(-j2*pi*tau_l*fs*k/K_bar) * aR_l * bT_l^T
%
%   SP positions are randomly generated in a 4λ x 4λ x 4λ cube.
%   Angles and ToA are derived from geometry (Eq. 49).
%
% Inputs:
%   params - system parameter struct
%
% Outputs:
%   Hk_all     - NR x NT x K array of channel matrices
%   alpha      - L x 1 complex gains ~ CN(0,1)
%   tau        - L x 1 time delays [s]
%   theta_az_R - L x 1 azimuth AoA [rad]
%   theta_el_R - L x 1 elevation AoA [rad]
%   theta_az_T - L x 1 azimuth AoD [rad]
%   theta_el_T - L x 1 elevation AoD [rad]
%   pl_all     - 3 x L SP positions [m]
%   d_c_R      - L x 1 distances SP->UT center [m]
%   d_c_T      - L x 1 distances SP->BS center [m]
%
% Runtime: O(L * (NR + NT) + NR*NT*K*L)
% =========================================================================

    L      = params.L;
    K      = params.K;
    K_bar  = params.K_bar;
    fs     = params.fs;
    lambda = params.lambda;
    c      = params.c;
    pT     = params.pT;
    pR     = params.pR;
    SP_range = params.SP_range;
    NR     = params.NR;
    NT     = params.NT;

    % -------------------------------------------------------
    % Generate L scattering point positions in near-field
    % A 4λ x 4λ x 4λ cube around the midpoint of BS-UT
    % -------------------------------------------------------
    midpoint = (pT + pR) / 2;
    pl_all = midpoint + SP_range * (rand(3, L) - 0.5);   % 3 x L

    % -------------------------------------------------------
    % Compute geometric parameters (Eq. 49)
    % -------------------------------------------------------
    theta_az_R = zeros(L,1);
    theta_el_R = zeros(L,1);
    theta_az_T = zeros(L,1);
    theta_el_T = zeros(L,1);
    tau        = zeros(L,1);
    d_c_R      = zeros(L,1);
    d_c_T      = zeros(L,1);

    for l = 1:L
        pl = pl_all(:,l);

        % Azimuth AoA: angle from UT to SP in xy-plane
        theta_az_R(l) = atan2(pl(2)-pR(2), pl(1)-pR(1)) + pi;

        % Elevation AoA
        dist_RL = norm(pl - pR);
        theta_el_R(l) = acos((pR(3)-pl(3)) / dist_RL);

        % Azimuth AoD
        theta_az_T(l) = atan2(pl(2)-pT(2), pl(1)-pT(1));

        % Elevation AoD
        dist_TL = norm(pl - pT);
        theta_el_T(l) = acos((pT(3)-pl(3)) / dist_TL);

        % ToA (sum of both legs)
        tau(l) = dist_RL/c + dist_TL/c;

        d_c_R(l) = dist_RL;
        d_c_T(l) = dist_TL;
    end

    % -------------------------------------------------------
    % Complex path gains ~ CN(0,1)
    % -------------------------------------------------------
    alpha = (randn(L,1) + 1j*randn(L,1)) / sqrt(2);

    % -------------------------------------------------------
    % Generate NR x NT channel matrices for each subcarrier
    % (Eq. 14)
    % -------------------------------------------------------
    Hk_all = zeros(NR, NT, K);

    % Precompute array response vectors for all L paths
    AR = zeros(NR, L);   % UT side
    BT = zeros(NT, L);   % BS side
    for l = 1:L
        [aR, bT] = near_field_array_response(params, pl_all(:,l));
        AR(:,l) = aR;
        BT(:,l) = bT;
    end

    % Subcarrier indices: k = 1, 2, ..., K  (selected, uniformly spaced)
    % Paper uses k as the index of selected subcarriers
    k_indices = round(linspace(1, K_bar, K));  % K selected from K_bar

    for ki = 1:K
        k = k_indices(ki);
        phase_delay = exp(-1j * 2*pi * tau .* fs * k / K_bar);  % L x 1
        % Hk = sum_l alpha_l * phase_l * aR_l * bT_l^T
        Hk = AR * diag(alpha .* phase_delay) * BT';   % NR x NT
        Hk_all(:,:,ki) = Hk;
    end
end
