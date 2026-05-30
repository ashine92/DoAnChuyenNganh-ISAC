function L_est = estimate_rank_mdl(Y, L_max)
% =========================================================================
% estimate_rank_mdl.m  (utils/)
% =========================================================================
% Description:
%   Estimates the CP rank (number of paths L) of a third-order tensor Y
%   using the Minimum Description Length (MDL) criterion.
%   Reference: Section 5.1 of paper, Eq. (26-27).
%
%   For each mode n, MDL minimizes:
%     MDL(l_n) = -2*J1*J2*J3/J_n * (J_n - l_n) *
%                 log[ prod_{zeta=l+1}^{J_n} eps_zeta^(1/(J_n-l_n))
%                      / (1/(J_n-l_n) * sum eps_zeta) ]
%                + l_n*(2*J_n - l_n) * log(J1*J2*J3/J_n)
%
%   Final estimate: L_hat = min over 3 modes of L_hat_n
%
% Inputs:
%   Y     - J1 x J2 x J3 tensor
%   L_max - maximum rank to search (default: min(J1,J2,J3)-1)
%
% Outputs:
%   L_est - estimated rank (number of paths)
%
% Runtime: O(J_n^2) per mode for eigendecomposition
% =========================================================================

    [J1, J2, J3] = size(Y);
    total_elem   = J1 * J2 * J3;

    if nargin < 2
        L_max = min([J1,J2,J3]) - 1;
    end

    L_hat_per_mode = zeros(1, 3);
    Jn_vec         = [J1, J2, J3];

    for n = 1:3
        Jn = Jn_vec(n);

        % Mode-n unfolding
        Y_n = tensor_unfold_mdl(Y, n);   % Jn x (total/Jn)

        % Sample covariance matrix eigenvalues
        Rn = (Y_n * Y_n') / (total_elem / Jn);
        eigenvalues = sort(real(eig(Rn)), 'descend');  % Jn x 1

        % MDL criterion (Eq. 27)
        mdl_vals = zeros(1, min(L_max, Jn-1));
        for ln = 0 : min(L_max, Jn-1) - 1
            tail_eigs = eigenvalues(ln+1:end);
            Ntail = length(tail_eigs);
            if Ntail <= 0 || any(tail_eigs <= 0)
                mdl_vals(ln+1) = inf;
                continue;
            end

            geom_mean  = prod(tail_eigs)^(1/Ntail);
            arith_mean = mean(tail_eigs);

            if geom_mean <= 0 || arith_mean <= 0
                mdl_vals(ln+1) = inf;
                continue;
            end

            log_ratio = log(geom_mean / arith_mean);  % <= 0
            pen_term  = ln * (2*Jn - ln) * log(total_elem / Jn);

            mdl_vals(ln+1) = -2 * (total_elem / Jn) * Ntail * log_ratio + pen_term;
        end

        [~, best_ln] = min(mdl_vals);
        L_hat_per_mode(n) = best_ln - 1;   % 0-indexed → actual rank
    end

    % Final estimate: minimum across modes (ensures uniqueness conditions)
    L_est = max(1, min(L_hat_per_mode));
    L_est = min(L_est, L_max);
end


%% --- Helper: tensor mode-n unfolding (local copy) -----------------------
function U = tensor_unfold_mdl(X, mode)
    sz    = size(X);
    ndims_X = length(sz);
    order = [mode, setdiff(1:ndims_X, mode)];
    X_perm = permute(X, order);
    U = reshape(X_perm, sz(mode), []);
end
