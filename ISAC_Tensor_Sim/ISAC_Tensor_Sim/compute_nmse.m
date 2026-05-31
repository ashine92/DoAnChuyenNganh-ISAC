function nmse_vals = compute_nmse(az_R, el_R, az_T, el_T, tau, pR, pl, ...
    az_R_h, el_R_h, az_T_h, el_T_h, tau_h, pR_h, pl_h)
% =========================================================================
% compute_nmse.m
% =========================================================================
% * Mathematical Background:
%   Computes the Normalized Mean Square Error (NMSE) for the channel 
%   parameters and localization.
%   NMSE(gamma) = ||gamma - gamma_hat||^2 / ||gamma||^2
%
% * Inputs:
%   az_R, el_R, az_T, el_T, tau, pR, pl - True parameters and locations
%   az_R_h, el_R_h, az_T_h, el_T_h, tau_h, pR_h, pl_h - Estimated values
%
% * Outputs:
%   nmse_vals - 1 x 7 array of NMSE values
%
% * MATLAB Implementation:
%   Applies an inline anonymous function for norm calculations.
%
% * Complexity Analysis:
%   O(L) where L is the number of paths/points.
% =========================================================================

    nmse = @(x, xh) norm(x(:)-xh(:))^2 / (norm(x(:))^2 + eps) * 100;

    nmse_vals = [ nmse(az_R, az_R_h), ...
                  nmse(el_R, el_R_h), ...
                  nmse(az_T, az_T_h), ...
                  nmse(el_T, el_T_h), ...
                  nmse(tau,  tau_h),  ...
                  nmse(pR,   pR_h),   ...
                  sum(sum((pl-pl_h).^2,1)) / (sum(sum(pl.^2,1)) + eps) * 100 ];
end
