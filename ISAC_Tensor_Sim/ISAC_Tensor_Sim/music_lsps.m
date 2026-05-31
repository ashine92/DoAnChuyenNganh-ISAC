function [az_R_hat, el_R_hat, az_T_hat, el_T_hat, tau_hat, pR_hat, pl_hat] = ...
    music_lsps(AR_hat, BT_hat, C_hat, params, pl_all, tau_true, az_R, el_R, az_T, el_T)
% =========================================================================
% music_lsps.m — MUSIC-Like Spectrum Peak Searching baseline
% =========================================================================
% Reference: [40] Pan et al., IEEE J. Sel. Top. SP 2023
%
% Algorithm:
%   1. Build covariance from AR_hat
%   2. 2D MUSIC spectrum search over (az, el) grid
%   3. Find L peaks → angle estimates
%   Grid resolution limits accuracy (inherent limitation vs proposed method)
% =========================================================================

    L      = params.L;
    NRy    = params.NRy;  NRz = params.NRz;
    lambda = params.lambda;
    d      = params.d;

    % --- AoA estimation via MUSIC ---
    R_hat = AR_hat * AR_hat';

    [U, D] = eig(R_hat);
    [~, sort_idx] = sort(real(diag(D)), 'descend');
    U_n = U(:, sort_idx(L+1:end));
    P_n = U_n * U_n';

    % Grid for 2D MUSIC search (coarser than proposed → inherent error)
    n_az = 90;  n_el = 45;
    az_grid = linspace(-pi, pi, n_az);
    el_grid = linspace(0.05, pi-0.05, n_el);

    [NY_R, NZ_R] = meshgrid(-(NRy-1)/2:(NRy-1)/2, -(NRz-1)/2:(NRz-1)/2);
    ny_R_vec = NY_R(:); nz_R_vec = NZ_R(:);

    music_spec = zeros(n_el, n_az);
    for ei = 1:n_el
        el = el_grid(ei);
        cos_el = cos(el);
        sin_el = sin(el);
        for ai = 1:n_az
            az = az_grid(ai);
            % Far-field steering vector (no near-field correction)
            phase = 2*pi*d/lambda * (ny_R_vec * sin(az)*sin_el + nz_R_vec * cos_el);
            a_ff = exp(-1j * phase);
            denom = real(a_ff' * P_n * a_ff);
            music_spec(ei,ai) = 1 / (denom + eps);
        end
    end

    % Find L peaks
    az_R_hat = zeros(L,1); el_R_hat = zeros(L,1);
    spec_copy = music_spec;
    for l = 1:L
        [~, idx] = max(spec_copy(:));
        [ei_max, ai_max] = ind2sub([n_el, n_az], idx);
        az_R_hat(l) = az_grid(ai_max);
        el_R_hat(l) = el_grid(ei_max);
        % Suppress neighborhood
        r = 3;
        ei_lo = max(1,ei_max-r); ei_hi = min(n_el,ei_max+r);
        ai_lo = max(1,ai_max-r); ai_hi = min(n_az,ai_max+r);
        spec_copy(ei_lo:ei_hi, ai_lo:ai_hi) = 0;
    end

    % Match to true paths using cost
    cost = abs(bsxfun(@minus, az_R_hat(:), az_R(:)'));
    perm = greedy_match(cost);
    az_R_hat = az_R_hat(perm);
    el_R_hat = el_R_hat(perm);

    % --- AoD estimation via MUSIC ---
    NTy = params.NTy; NTz = params.NTz;
    [NY_T, NZ_T] = meshgrid(-(NTy-1)/2:(NTy-1)/2, -(NTz-1)/2:(NTz-1)/2);
    ny_T_vec = NY_T(:); nz_T_vec = NZ_T(:);

    R_T = BT_hat * BT_hat';
    [U_T, D_T] = eig(R_T);
    [~, si_T] = sort(real(diag(D_T)), 'descend');
    U_n_T = U_T(:, si_T(L+1:end));
    P_n_T = U_n_T * U_n_T';

    spec_T = zeros(n_el, n_az);
    for ei = 1:n_el
        el = el_grid(ei);
        cos_el = cos(el); sin_el = sin(el);
        for ai = 1:n_az
            az = az_grid(ai);
            phase = 2*pi*d/lambda * (ny_T_vec*sin(az)*sin_el + nz_T_vec*cos_el);
            a_ff = exp(-1j*phase);
            denom = real(a_ff' * P_n_T * a_ff);
            spec_T(ei,ai) = 1 / (denom + eps);
        end
    end

    az_T_hat = zeros(L,1); el_T_hat = zeros(L,1);
    for l = 1:L
        [~,idx] = max(spec_T(:));
        [ei_max,ai_max] = ind2sub([n_el,n_az],idx);
        az_T_hat(l) = az_grid(ai_max); el_T_hat(l) = el_grid(ei_max);
        r=3;
        spec_T(max(1,ei_max-r):min(n_el,ei_max+r), max(1,ai_max-r):min(n_az,ai_max+r)) = 0;
    end

    cost_T = abs(bsxfun(@minus, az_T_hat(:), az_T(:)'));
    perm_T = greedy_match(cost_T);
    az_T_hat = az_T_hat(perm_T);
    el_T_hat = el_T_hat(perm_T);

    % --- ToA ---
    tau_hat = estimate_toa(C_hat, params);
    cost_tau = abs(bsxfun(@minus, tau_hat(:), tau_true(:)'));
    perm_tau = greedy_match(cost_tau);
    tau_hat = tau_hat(perm_tau);

    % --- Localization ---
    pR_hat = localize_ut(tau_hat, az_R_hat, el_R_hat, az_T_hat, el_T_hat, params);
    pl_hat = localize_sps(tau_hat, az_R_hat, el_R_hat, az_T_hat, el_T_hat, pR_hat, params);
end


function perm = greedy_match(cost)
    L = size(cost, 1);
    perm = zeros(1, L);
    used = false(1, L);
    for j = 1:L
        c = cost(:,j);
        c(used) = inf;
        [~, i] = min(c);
        perm(j) = i;
        used(i) = true;
    end
end
