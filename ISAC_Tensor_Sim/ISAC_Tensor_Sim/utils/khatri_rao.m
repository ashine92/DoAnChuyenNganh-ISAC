function KR = khatri_rao(A, B)
% =========================================================================
% khatri_rao.m  (utils/)
% =========================================================================
% Description:
%   Computes the Khatri-Rao (column-wise Kronecker) product of two matrices.
%
%   For A ∈ C^{m x L} and B ∈ C^{n x L}:
%     A ⊙ B ∈ C^{(m*n) x L}
%     (A ⊙ B)_l = kron(A_l, B_l)  for each column l
%
% Inputs:
%   A - m x L matrix
%   B - n x L matrix
%
% Output:
%   KR - (m*n) x L Khatri-Rao product
%
% Runtime: O(m*n*L)
% =========================================================================

    [m, L] = size(A);
    n      = size(B, 1);

    KR = zeros(m*n, L, class(A));

    for l = 1:L
        KR(:,l) = kron(A(:,l), B(:,l));
    end
end
