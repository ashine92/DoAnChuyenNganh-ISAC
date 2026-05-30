function [A_hat, B_hat, C_hat, iter_count, fit_history] = cp_als(Y, L, params)
% =========================================================================
% cp_als.m  —  CP Decomposition via Alternating Least Squares
% =========================================================================
% Model: Y ≈ Σ_l a_l ∘ b_l ∘ c_l   (Eq. 28)
%
% ALS subproblem for each factor (Kolda & Bader 2009, §4.2):
%   min_A ||Y_(1) - A*(C⊙B)^T||^2
%   Solution:  A^T = pinv(C⊙B) * Y_(1)^T   [minimum-norm LS]
%
% Using pinv(C⊙B) instead of \ or Gram-matrix inversion is numerically
% stable even when C⊙B is rank-deficient (small F,T,K or near-collinear paths).
%
% Tensor unfoldings (MATLAB column-major, matches paper Eq. 25):
%   Y_(1) = reshape(permute(Y,[1,2,3]), J1, J2*J3)   j2 varies fast
%   Y_(2) = reshape(permute(Y,[2,1,3]), J2, J1*J3)   j1 varies fast
%   Y_(3) = reshape(permute(Y,[3,1,2]), J3, J1*J2)   j1 varies fast
%
% Reconstruction (inverse of mode-3 unfold):
%   Y3_hat = C*(B⊙A)^T  →  permute(reshape(Y3_hat,[J3,J1,J2]),[2,3,1])
% =========================================================================

    max_iter = params.max_iter;
    eps_tol  = params.eps_tol;

    [J1, J2, J3] = size(Y);
    norm_Y = max(eps, norm(Y(:)));

    % --- Mode unfoldings ---
    Y1 = reshape(permute(Y,[1,2,3]), J1, J2*J3);   % J1 x (J2*J3)
    Y2 = reshape(permute(Y,[2,1,3]), J2, J1*J3);   % J2 x (J1*J3)
    Y3 = reshape(permute(Y,[3,1,2]), J3, J1*J2);   % J3 x (J1*J2)

    % --- Initialise: leading L left singular vectors (scaled by sqrt(σ)) ---
    A_hat = svd_init_als(Y1, L);
    B_hat = svd_init_als(Y2, L);
    C_hat = svd_init_als(Y3, L);

    fit_history = zeros(max_iter, 1);
    iter_count  = 0;
    prev_fit    = inf;

    % --- ALS iterations ---
    for m = 1:max_iter

        % Update A:  A^T = pinv(C⊙B) * Y1^T
        CB    = kr_als(C_hat, B_hat);          % (J3*J2) x L
        A_hat = (pinv(CB) * Y1')';             % J1 x L

        % Update B:  B^T = pinv(C⊙A) * Y2^T
        CA    = kr_als(C_hat, A_hat);          % (J3*J1) x L
        B_hat = (pinv(CA) * Y2')';             % J2 x L

        % Update C:  C^T = pinv(B⊙A) * Y3^T
        BA    = kr_als(B_hat, A_hat);          % (J2*J1) x L
        C_hat = (pinv(BA) * Y3')';             % J3 x L

        % Fit
        Y_hat    = cp_reconstruct_als(A_hat, B_hat, C_hat, J1, J2, J3);
        residual = norm(Y(:) - Y_hat(:)) / norm_Y;
        fit_history(m) = residual;
        iter_count = m;

        if m > 5 && abs(prev_fit - residual) < eps_tol
            break;
        end
        prev_fit = residual;
    end

    fit_history = fit_history(1:iter_count);
end


function KR = kr_als(A, B)
% Khatri-Rao product: KR(:,l) = kron(A(:,l), B(:,l))
% A: m×L,  B: n×L  →  KR: (m*n)×L  — fully vectorised
    [m, L] = size(A);
    n  = size(B, 1);
    % Expand and multiply: A(i,l)*B(j,l) at row i+m*(j-1)
    KR = reshape( ...
            bsxfun(@times, reshape(A,[m,1,L]), reshape(B,[1,n,L])), ...
            m*n, L);
end


function U = svd_init_als(M, L)
% First L left singular vectors, scaled by sqrt(singular values).
    k = max(1, min(L, min(size(M))-1));
    [U, S, ~] = svds(M, k);
    s = max(eps, sqrt(diag(S)))';   % 1×k scale
    U = U .* s;
    if size(U,2) < L
        pad = (randn(size(U,1), L-size(U,2)) + ...
               1j*randn(size(U,1), L-size(U,2))) * norm(U,'fro') / (L*size(U,1));
        U = [U, pad];
    end
end


function Y_hat = cp_reconstruct_als(A, B, C, J1, J2, J3)
% Y_hat = Σ_l a_l ∘ b_l ∘ c_l  via inverse mode-3 unfolding.
%   Y_(3) = reshape(permute(Y,[3,1,2]), J3, J1*J2)   j1 fast, j2 slow
%   Y_(3) = C*(B⊙A)^T  with B⊙A = kron(b_l,a_l)  j1 fast ✓
%   Inverse: reshape Y3_hat to [J3,J1,J2] (col-major), then permute([2,3,1])
    BA    = kr_als(B, A);                           % (J2*J1)×L  j1 fast
    Y3hat = C * BA';                                 % J3×(J2*J1)
    Y_hat = permute(reshape(Y3hat, [J3,J1,J2]), [2,3,1]);
end
