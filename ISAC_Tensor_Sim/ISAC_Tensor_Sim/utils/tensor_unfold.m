function U = tensor_unfold(X, mode)
% =========================================================================
% tensor_unfold.m  (utils/)
% =========================================================================
% Description:
%   Computes the mode-n matricization (unfolding/flattening) of a
%   third-order tensor X ∈ C^{J1 x J2 x J3}.
%
%   The mode-n unfolding [X]_(n) ∈ C^{J_n x (J_{n+1}*...*J_N*J_1*...*J_{n-1})}
%   arranges the mode-n fibers as columns.
%
%   Ordering convention (consistent with Khatri-Rao in cp_als.m):
%     mode=1: U ∈ C^{J1 x (J2*J3)}, fibers ordered j2 fast, j3 slow
%     mode=2: U ∈ C^{J2 x (J1*J3)}, fibers ordered j1 fast, j3 slow
%     mode=3: U ∈ C^{J3 x (J1*J2)}, fibers ordered j1 fast, j2 slow
%
%   This matches the unfolding convention used in Eq. (5) and (25) of the paper:
%     [Y]_(1) = A * (C ⊙ B)^T   → columns indexed by (j2, j3) with j2 fast
%
% Inputs:
%   X    - J1 x J2 x J3 complex tensor
%   mode - unfolding mode (1, 2, or 3)
%
% Outputs:
%   U    - J_n x (J1*J2*J3/J_n) unfolded matrix
%
% Runtime: O(J1*J2*J3) (permute + reshape)
% =========================================================================

    sz      = size(X);
    ndims_X = ndims(X);

    % Move target mode to first position, then reshape
    order   = [mode, setdiff(1:ndims_X, mode)];
    X_perm  = permute(X, order);
    U       = reshape(X_perm, sz(mode), []);
end
