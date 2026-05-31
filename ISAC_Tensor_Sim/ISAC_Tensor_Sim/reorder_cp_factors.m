function [A_reordered, B_reordered, C_reordered, reorder_idx] = reorder_cp_factors(A_hat, B_hat, C_hat, params)
% =========================================================================
% reorder_cp_factors.m - Fix CP decomposition path ordering ambiguity
% =========================================================================
% 
% PROBLEM: CP-ALS converges to a valid solution, but with paths permuted
% because the tensor model is permutation-invariant in the factor indices.
%
% SOLUTION: Reorder factors based on canonical criterion (sorted ToA)
% Uses robust 1D correlation search via estimate_toa.
% =========================================================================

    % Estimate ToA using the robust 1D correlation search
    tau_hat = estimate_toa(C_hat, params);
    
    % Sort by estimated ToA to get canonical ordering
    [~, reorder_idx] = sort(tau_hat);
    
    % Reorder factors
    A_reordered = A_hat(:, reorder_idx);
    B_reordered = B_hat(:, reorder_idx);
    C_reordered = C_hat(:, reorder_idx);
end
