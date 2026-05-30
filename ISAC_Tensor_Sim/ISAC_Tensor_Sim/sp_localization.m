function pl_hat = sp_localization(tau_hat, theta_az_R_hat, theta_el_R_hat, ...
                                   theta_az_T_hat, theta_el_T_hat, pR_hat, params)
% =========================================================================
% sp_localization.m
% =========================================================================
% Description:
%   Estimates the 3D positions of Scattering Points (SPs) from estimated
%   channel parameters and the estimated UT position.
%
%   Mathematical model (Eq. 46-47):
%
%   Distance from SP to BS along gT,l:
%     d_c_T,l = gT,l^T * (pl - pT) / ||gT,l||^2
%
%   Distance from SP to UT along gR,l:
%     d_c_R,l = gR,l^T * (pl - pR_hat) / ||gR,l||^2
%
%   SP position from intersection (Eq. 47):
%     pl_hat = (Q_T,l + Q_R,l)^{-1} * (Q_T,l * pT + Q_R,l * pR_hat)
%   where:
%     Q_T,l = I3 - gT,l * gT,l^T
%     Q_R,l = I3 - gR,l * gR,l^T
%
% Inputs:
%   tau_hat         - L x 1 estimated ToAs [s]
%   theta_az_R_hat  - L x 1 estimated azimuth AoA [rad]
%   theta_el_R_hat  - L x 1 estimated elevation AoA [rad]
%   theta_az_T_hat  - L x 1 estimated azimuth AoD [rad]
%   theta_el_T_hat  - L x 1 estimated elevation AoD [rad]
%   pR_hat          - 3 x 1 estimated UT position [m]
%   params          - struct with fields: pT, L
%
% Outputs:
%   pl_hat - 3 x L estimated SP positions [m]
%
% Runtime: O(L * 27)  (matrix operations)
% =========================================================================

    L   = params.L;
    pT  = params.pT;
    I3  = eye(3);

    pl_hat = zeros(3, L);

    for l = 1:L
        az_T = theta_az_T_hat(l);
        el_T = theta_el_T_hat(l);
        az_R = theta_az_R_hat(l);
        el_R = theta_el_R_hat(l);

        % Direction vectors (unit, 3D)
        gT = [cos(az_T)*sin(el_T); sin(az_T)*sin(el_T); cos(el_T)];
        gR = [cos(az_R)*sin(el_R); sin(az_R)*sin(el_R); cos(el_R)];

        gT = gT / (norm(gT) + eps);
        gR = gR / (norm(gR) + eps);

        % Projection matrices (Eq. 47)
        Q_T = I3 - gT * gT';
        Q_R = I3 - gR * gR';

        % SP position (Eq. 47)
        A = Q_T + Q_R;
        b = Q_T * pT + Q_R * pR_hat;

        pl_hat(:,l) = A \ b;
    end
end
