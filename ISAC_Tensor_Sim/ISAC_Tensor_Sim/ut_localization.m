function pR_hat = ut_localization(tau_hat, theta_az_R_hat, theta_el_R_hat, ...
                                   theta_az_T_hat, theta_el_T_hat, params)
% =========================================================================
% ut_localization.m
% =========================================================================
% Description:
%   Estimates the 3D position of the User Terminal (UT) from estimated
%   channel parameters using geometric constraints.
%
%   Mathematical model (Eq. 42-45):
%     For each path l, define:
%       gT,l = [cos(az_T)*sin(el_T), sin(az_T)*sin(el_T), cos(el_T)]^T  (tx dir)
%       gR,l = [cos(az_R)*sin(el_R), sin(az_R)*sin(el_R), cos(el_R)]^T  (rx dir)
%
%     Line equation: pR = pT + c*tau_l*upsilon_l*gT,l - c*tau_l*(1-upsilon_l)*gR,l
%                       = eta_l + upsilon_l * u_l
%     where eta_l = pT - c*tau_l*gR,l, u_l = c*tau_l*(gT,l + gR,l)
%
%   UT position minimizes (Eq. 44):
%     L(pR) = sum_l xi_l * ||pR - (eta_l + u_l^T(pR-eta_l)*u_l)||^2_2
%
%   Closed-form solution (Eq. 45):
%     pR_hat = [sum_l xi_l*(I3 - u_bar_l*u_bar_l^T)]^{-1}
%               * sum_l xi_l*(I3 - u_bar_l*u_bar_l^T)*eta_l
%
% Inputs:
%   tau_hat         - L x 1 estimated ToAs [s]
%   theta_az_R_hat  - L x 1 estimated azimuth AoA [rad]
%   theta_el_R_hat  - L x 1 estimated elevation AoA [rad]
%   theta_az_T_hat  - L x 1 estimated azimuth AoD [rad]
%   theta_el_T_hat  - L x 1 estimated elevation AoD [rad]
%   params          - struct with fields: pT, c, L
%
% Outputs:
%   pR_hat - 3 x 1 estimated UT position [m]
%
% Runtime: O(L * 9)  (matrix operations on 3x3 matrices)
% =========================================================================

    L   = params.L;
    pT  = params.pT;
    c   = params.c;

    I3  = eye(3);
    A_sum = zeros(3, 3);
    b_sum = zeros(3, 1);

    % Equal weights xi_l = 1 for all paths
    xi = ones(L, 1);

    for l = 1:L
        az_T = theta_az_T_hat(l);
        el_T = theta_el_T_hat(l);
        az_R = theta_az_R_hat(l);
        el_R = theta_el_R_hat(l);

        % Direction vectors (unit, 3D)
        gT = [cos(az_T)*sin(el_T); sin(az_T)*sin(el_T); cos(el_T)];
        gR = [cos(az_R)*sin(el_R); sin(az_R)*sin(el_R); cos(el_R)];

        % Normalize just in case
        gT = gT / (norm(gT) + eps);
        gR = gR / (norm(gR) + eps);

        % eta_l = pT - c*tau_l*gR (Eq. 42)
        eta_l = pT - c * tau_hat(l) * gR;

        % u_l = c*tau_l*(gT + gR)  (Eq. 42)
        u_l = c * tau_hat(l) * (gT + gR);

        % Normalized u_bar_l = u_l / ||u_l||
        u_bar = u_l / (norm(u_l) + eps);

        % Projection matrix: I3 - u_bar * u_bar^T
        P_l = I3 - u_bar * u_bar';

        % Accumulate (Eq. 45)
        A_sum = A_sum + xi(l) * P_l;
        b_sum = b_sum + xi(l) * P_l * eta_l;
    end

    % Solve linear system: A_sum * pR = b_sum
    pR_hat = A_sum \ b_sum;
end
