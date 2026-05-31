function pl_hat = localize_sps(tau_hat, theta_az_R_hat, theta_el_R_hat, ...
                                   theta_az_T_hat, theta_el_T_hat, pR_hat, params)
% =========================================================================
% localize_sps.m
% =========================================================================
% Estimates SP positions using line intersection (Eq. 46-47).
% Uses the cos(az) sign convention: az ∈ [-π/2, π/2] from ESPRIT,
% resolves signs using the constraint that SPs lie between BS and UT.
% =========================================================================

    L   = params.L;
    pT  = params.pT;
    c   = params.c;
    I3  = eye(3);

    pl_hat = zeros(3, L);

    for l = 1:L
        sin_az_T = sin(theta_az_T_hat(l));
        sin_az_R = sin(theta_az_R_hat(l));
        sin_el_T = sin(theta_el_T_hat(l));
        sin_el_R = sin(theta_el_R_hat(l));
        cos_el_T = cos(theta_el_T_hat(l));
        cos_el_R = cos(theta_el_R_hat(l));

        % Try both signs for cos(az_R) and keep cos(az_T) > 0
        best_err = inf;
        best_pl = zeros(3,1);

        for sign_R = [-1, 1]
            cos_az_R = sign_R * cos(theta_az_R_hat(l));
            cos_az_T = cos(theta_az_T_hat(l));  % positive

            gT = [cos_az_T*sin_el_T; sin_az_T*sin_el_T; cos_el_T];
            gR = [cos_az_R*sin_el_R; sin_az_R*sin_el_R; cos_el_R];
            gT = gT / (norm(gT)+eps);
            gR = gR / (norm(gR)+eps);

            Q_T = I3 - gT * gT';
            Q_R = I3 - gR * gR';

            A = Q_T + Q_R + eye(3)*1e-10;
            b = Q_T * pT + Q_R * pR_hat;
            pl_try = A \ b;

            % Check: SP should be reachable from both arrays
            d_T = norm(pl_try - pT);
            d_R = norm(pl_try - pR_hat);
            tau_check = (d_T + d_R) / c;
            err = abs(tau_check - tau_hat(l));

            if err < best_err
                best_err = err;
                best_pl = pl_try;
            end
        end

        pl_hat(:,l) = best_pl;
    end
end
