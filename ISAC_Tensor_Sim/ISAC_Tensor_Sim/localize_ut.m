function pR_hat = localize_ut(tau_hat, theta_az_R_hat, theta_el_R_hat, ...
                                   theta_az_T_hat, theta_el_T_hat, params)
% =========================================================================
% localize_ut.m
% =========================================================================
% Estimates the 3D position of the User Terminal (UT).
%
% Input angles are in [-π/2, π/2] for azimuth (cos(az) ≥ 0).
% We try both signs of cos(az) for each path and pick the combination
% that gives the most consistent localization.
%
% Method: Eq. 42-45 from paper
%   pR = [Σ ξ_l (I - ū_l ū_l^T)]^(-1) * [Σ ξ_l (I - ū_l ū_l^T) η_l]
% =========================================================================

    L   = params.L;
    pT  = params.pT;
    c   = params.c;
    I3  = eye(3);

    % For each path, resolve the cos(az) sign ambiguity
    % by trying both signs and selecting the one where the SP
    % lies roughly between pT and pR (estimated iteratively)

    % Step 1: Get initial pR estimate using cos(az) > 0 for all paths
    pR_best = solve_localization(tau_hat, theta_az_R_hat, theta_el_R_hat, ...
                                  theta_az_T_hat, theta_el_T_hat, pT, c, L, I3, ...
                                  ones(L,1), ones(L,1));

    % Step 2: Try flipping cos(az) signs for each path
    % For AoD (from BS), the SP should have x > pT_x typically,
    % so cos(az_T) > 0 is usually correct.
    % For AoA (from UT), try sign = -1 (cos(az_R) < 0, sp is behind UT in x)
    % This is more likely given the geometry (pR_x > midpoint_x)
    best_residual = inf;
    best_signs_R = ones(L,1);
    best_signs_T = ones(L,1);

    % Try all 2^L sign combinations for AoA (keep AoD positive)
    for mask = 0:(2^L - 1)
        signs_R = ones(L, 1);
        for l = 1:L
            if bitget(mask, l)
                signs_R(l) = -1;
            end
        end
        signs_T = ones(L, 1);

        pR_try = solve_localization(tau_hat, theta_az_R_hat, theta_el_R_hat, ...
                                     theta_az_T_hat, theta_el_T_hat, pT, c, L, I3, ...
                                     signs_R, signs_T);

        % Residual: check consistency of all paths
        res = compute_residual(pR_try, pT, tau_hat, theta_az_R_hat, theta_el_R_hat, ...
                               theta_az_T_hat, theta_el_T_hat, c, L, signs_R, signs_T);

        if res < best_residual
            best_residual = res;
            pR_best = pR_try;
            best_signs_R = signs_R;
            best_signs_T = signs_T;
        end
    end

    pR_hat = pR_best;
end


function pR_est = solve_localization(tau, az_R, el_R, az_T, el_T, ...
                                      pT, c, L, I3, signs_R, signs_T)
    A_sum = zeros(3,3);
    b_sum = zeros(3,1);
    xi = ones(L,1);

    for l = 1:L
        sin_az_R = sin(az_R(l));
        cos_az_R = signs_R(l) * cos(az_R(l));  % Apply sign
        sin_el_R = sin(el_R(l));
        cos_el_R = cos(el_R(l));

        sin_az_T = sin(az_T(l));
        cos_az_T = signs_T(l) * cos(az_T(l));
        sin_el_T = sin(el_T(l));
        cos_el_T = cos(el_T(l));

        gT = [cos_az_T*sin_el_T; sin_az_T*sin_el_T; cos_el_T];
        gR = [cos_az_R*sin_el_R; sin_az_R*sin_el_R; cos_el_R];

        gT = gT / (norm(gT) + eps);
        gR = gR / (norm(gR) + eps);

        eta_l = pT - c * tau(l) * gR;
        u_l = c * tau(l) * (gT + gR);

        u_norm = norm(u_l);
        if u_norm < eps, continue; end
        u_bar = u_l / u_norm;

        P_l = I3 - u_bar * u_bar';

        A_sum = A_sum + xi(l) * P_l;
        b_sum = b_sum + xi(l) * P_l * eta_l;
    end

    A_sum = A_sum + eye(3) * 1e-10;
    pR_est = A_sum \ b_sum;
end


function res = compute_residual(pR_est, pT, tau, az_R, el_R, az_T, el_T, ...
                                 c, L, signs_R, signs_T)
    res = 0;
    for l = 1:L
        sin_az_R = sin(az_R(l));
        cos_az_R = signs_R(l) * cos(az_R(l));
        sin_az_T = sin(az_T(l));
        cos_az_T = signs_T(l) * cos(az_T(l));

        gT = [cos_az_T*sin(el_T(l)); sin_az_T*sin(el_T(l)); cos(el_T(l))];
        gR = [cos_az_R*sin(el_R(l)); sin_az_R*sin(el_R(l)); cos(el_R(l))];
        gT = gT / (norm(gT)+eps);
        gR = gR / (norm(gR)+eps);

        % SP from BS side: pl_T = pT + d_T * gT
        % SP from UT side: pl_R = pR + d_R * gR
        % Should match with d_T + d_R = c*tau
        % Residual: ||pl_T - pl_R||² for best d_T, d_R

        % Least squares: pT + d_T*gT = pR + d_R*gR
        % => [gT, -gR] * [d_T; d_R] = pR_est - pT
        A = [gT, -gR];
        b_vec = pR_est - pT;
        x = pinv(A) * b_vec;  % [d_T; d_R]
        d_T = x(1); d_R = x(2);

        % Residual: consistency check
        pl_T = pT + d_T * gT;
        pl_R = pR_est + d_R * gR;
        res = res + norm(pl_T - pl_R)^2;
        % Also penalize negative distances or wrong total
        if d_T < 0 || d_R < 0
            res = res + 100;
        end
        res = res + (d_T + d_R - c*tau(l))^2;
    end
end
