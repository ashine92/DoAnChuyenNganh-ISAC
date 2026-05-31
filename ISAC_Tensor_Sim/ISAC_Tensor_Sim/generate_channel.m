function [Hk_all, alpha, tau, theta_az_R, theta_el_R, theta_az_T, theta_el_T, pl_all, d_c_R, d_c_T] = ...
    generate_channel(params)
% =========================================================================
% generate_channel.m
% =========================================================================
% Generates the near-field mmWave MIMO-OFDM channel for all K selected
% subcarriers using NLoS-only model with spherical wavefront.
%
% Frequency-domain channel at subcarrier k (Eq. 14):
%   Hk = sum_l alpha_l * exp(-j2*pi*tau_l*fs*k/K_bar) * aR_l * bT_l^T
%
% Angle convention (consistent with paper Eq. 49):
%   AoA: direction from UT to SP
%     azimuth: atan2(pl_y - pR_y, pl_x - pR_x)
%     elevation: acos((pl_z - pR_z) / d_R)  (angle from z-axis)
%   AoD: direction from BS to SP
%     azimuth: atan2(pl_y - pT_y, pl_x - pT_x)
%     elevation: acos((pl_z - pT_z) / d_T)  (angle from z-axis)
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
    % -------------------------------------------------------
    midpoint = (pT + pR) / 2;
    pl_all = midpoint + SP_range .* (rand(3, L) - 0.5);   % 3 x L
    
    % FIX: A planar array in the y-z plane cannot resolve +x vs -x.
    % Force the x-coordinate to be strictly positive to match the 
    % assumed azimuth range [-pi/2, pi/2] (i.e., in front of the array).
    pl_all(1, :) = abs(pl_all(1, :)) + 2; % ensure x >= 2
    
    % -------------------------------------------------------
    % Compute geometric parameters
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

        % Distances
        dist_RL = norm(pl - pR);
        dist_TL = norm(pl - pT);

        % Elevation: angle from z-axis
        theta_el_R(l) = acos(max(-1, min(1, (pl(3)-pR(3)) / dist_RL)));
        theta_el_T(l) = acos(max(-1, min(1, (pl(3)-pT(3)) / dist_TL)));

        % Observable quantities for yoz-plane UPA:
        %   β = sin(az)*sin(el) = Δy / dist
        %   σ = cos(el) = Δz / dist
        % Since UPA cannot resolve sign of cos(az),
        % define az = atan2(sin(az), |cos(az)|) → az ∈ [-π/2, π/2]
        sin_el_R = sin(theta_el_R(l));
        if abs(sin_el_R) > 1e-10
            sin_az_R = (pl(2)-pR(2)) / (dist_RL * sin_el_R);
            sin_az_R = max(-1, min(1, sin_az_R));
            cos_az_R = sqrt(max(0, 1 - sin_az_R^2));  % always non-negative
            theta_az_R(l) = atan2(sin_az_R, cos_az_R);
        else
            theta_az_R(l) = 0;
        end

        sin_el_T = sin(theta_el_T(l));
        if abs(sin_el_T) > 1e-10
            sin_az_T = (pl(2)-pT(2)) / (dist_TL * sin_el_T);
            sin_az_T = max(-1, min(1, sin_az_T));
            cos_az_T = sqrt(max(0, 1 - sin_az_T^2));
            theta_az_T(l) = atan2(sin_az_T, cos_az_T);
        else
            theta_az_T(l) = 0;
        end

        % ToA (sum of both legs)
        tau(l) = (dist_RL + dist_TL) / c;

        d_c_R(l) = dist_RL;
        d_c_T(l) = dist_TL;
    end

    % -------------------------------------------------------
    % Complex path gains ~ CN(0,1)
    % -------------------------------------------------------
    alpha = (randn(L,1) + 1j*randn(L,1)) / sqrt(2);

    % -------------------------------------------------------
    % Generate NR x NT channel matrices for each subcarrier
    % -------------------------------------------------------
    Hk_all = zeros(NR, NT, K);

    % Precompute array response vectors for all L paths
    AR = zeros(NR, L);   % UT side
    BT = zeros(NT, L);   % BS side
    for l = 1:L
        [aR, bT] = nearfield_array_response(params, pl_all(:,l));
        AR(:,l) = aR;
        BT(:,l) = bT;
    end

    % Subcarrier indices: K selected uniformly from K_bar
    k_indices = round(linspace(1, K_bar, K));

    for ki = 1:K
        k = k_indices(ki);
        % phase_delay corresponds to delay tau at subcarrier index k
        phase_delay = exp(-1j * 2*pi * tau * params.fs .* k / K_bar);  % L x 1
        Hk = AR * diag(alpha .* phase_delay) * BT';   % NR x NT
        Hk_all(:,:,ki) = Hk;
    end
end
