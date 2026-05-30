function nmse = compute_nmse(x_true, x_hat)
% =========================================================================
% compute_nmse.m  (utils/)
% =========================================================================
% Description:
%   Computes the Normalized Mean Square Error (NMSE) between true and
%   estimated parameter vectors.
%
%   NMSE(γ) = ||γ - γ_hat||²₂ / ||γ||²₂
%
%   For SP localization (multiple targets):
%   NMSE(p_l) = Σ_l ||p_l - p_hat_l||²₂ / Σ_l ||p_l||²₂
%
% Inputs:
%   x_true - true parameter vector / matrix (real or complex)
%   x_hat  - estimated parameter vector / matrix
%
% Outputs:
%   nmse   - scalar NMSE value (non-negative)
%
% Usage example:
%   nmse_aoa = compute_nmse(theta_az_R_true, theta_az_R_hat);
%   nmse_sp  = compute_nmse(pl_all, pl_hat);   % 3 x L matrices
% =========================================================================

    x_true = x_true(:);
    x_hat  = x_hat(:);

    norm_true = norm(x_true);
    if norm_true < eps
        nmse = norm(x_true - x_hat)^2;   % Unnormalized when true is zero
    else
        nmse = norm(x_true - x_hat)^2 / norm_true^2;
    end
end
