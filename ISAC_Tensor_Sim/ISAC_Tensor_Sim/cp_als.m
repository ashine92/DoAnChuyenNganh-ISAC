function [A_hat, B_hat, C_hat, iter_count, fit_history] = cp_als(Y, L, params)
% =========================================================================
% cp_als.m  —  CP Decomposition via Alternating Least Squares
% =========================================================================
% Model: Y ≈ Σ_l a_l ∘ b_l ∘ c_l   (Eq. 28)
%
% ALS update (Kolda & Bader 2009, §4.2):
%   Mode-1: Y_(1) = A * (C⊙B)^T   =>  A^T = pinv(C⊙B) * Y_(1)^T
%   Mode-2: Y_(2) = B * (C⊙A)^T   =>  B^T = pinv(C⊙A) * Y_(2)^T
%   Mode-3: Y_(3) = C * (B⊙A)^T   =>  C^T = pinv(B⊙A) * Y_(3)^T
%
% Uses pinv for numerical stability when Gram matrix is ill-conditioned
% (which happens when F > NR or T > NT).
% =========================================================================

    max_iter = params.max_iter;
    eps_tol  = params.eps_tol;

    [J1, J2, J3] = size(Y);
    norm_Y = max(eps, norm(Y(:)));

    % --- Mode unfoldings ---
    Y1 = reshape(Y, J1, J2*J3);                      % J1 x (J2*J3)
    Y2 = reshape(permute(Y, [2,1,3]), J2, J1*J3);    % J2 x (J1*J3)
    Y3 = reshape(permute(Y, [3,1,2]), J3, J1*J2);    % J3 x (J1*J2)

    % --- Initialise: Support multiple methods ---
    init_method = 'svd';  % Default: SVD-based
    if isfield(params, 'cp_als_init_method')
        init_method = params.cp_als_init_method;
    end

    switch lower(init_method)
        case 'random'
            % Random Gaussian initialization
            A_hat = (randn(J1, L) + 1j*randn(J1, L)) / sqrt(2);
            B_hat = (randn(J2, L) + 1j*randn(J2, L)) / sqrt(2);
            C_hat = (randn(J3, L) + 1j*randn(J3, L)) / sqrt(2);

        case 'svd'
            % SVD-based initialization (structured)
            A_hat = svd_init(Y1, L);
            B_hat = svd_init(Y2, L);
            C_hat = svd_init(Y3, L);

        case 'data_driven'
            % Data-driven: use strongest SVD directions + noise
            A_hat = svd_init(Y1, L);
            B_hat = svd_init(Y2, L);
            C_hat = svd_init(Y3, L);
            % Add small data-driven boost from unfoldings
            [U, ~, ~] = svds(Y2, L);  % Extra insight from mode-2
            A_hat = 0.7*A_hat + 0.3*U;

        case 'ensemble'
            % Will be handled outside this function
            % Placeholder: use SVD initially
            A_hat = svd_init(Y1, L);
            B_hat = svd_init(Y2, L);
            C_hat = svd_init(Y3, L);

        otherwise
            error('Unknown CP-ALS initialization method: %s', init_method);
    end

    fit_history = zeros(max_iter, 1);
    iter_count  = 0;
    prev_err    = inf;

    % --- ALS iterations ---
    for m = 1:max_iter

        % Update A:  Y_(1) = A * (C⊙B)^T  =>  A = Y_(1) * pinv((C⊙B)^T)
        CB    = kr(C_hat, B_hat);             % (J3*J2) x L
        A_hat = (pinv(CB) * Y1.').';          % J1 x L

        % Update B:  Y_(2) = B * (C⊙A)^T
        CA    = kr(C_hat, A_hat);             % (J3*J1) x L
        B_hat = (pinv(CA) * Y2.').';          % J2 x L

        % Update C:  Y_(3) = C * (B⊙A)^T
        BA    = kr(B_hat, A_hat);             % (J2*J1) x L
        C_hat = (pinv(BA) * Y3.').';          % J3 x L

        % Column-wise normalization for numerical stability
        % Only normalize to prevent overflow/underflow, don't collapse all magnitude to C
        for l = 1:L
            nA = norm(A_hat(:,l));
            nB = norm(B_hat(:,l));
            nC = norm(C_hat(:,l));
            
            % Geometric normalization: distribute magnitude across factors
            if nA > eps && nB > eps && nC > eps
                % Normalize each factor by cube root of product (keeps magnitude)
                overall_scale = (nA * nB * nC)^(1/3);
                
                A_hat(:,l) = A_hat(:,l) / (nA / overall_scale);
                B_hat(:,l) = B_hat(:,l) / (nB / overall_scale);
                C_hat(:,l) = C_hat(:,l) / (nC / overall_scale);
            end
        end

        % Compute fit (relative error)
        Y_hat = reconstruct(A_hat, B_hat, C_hat, J1, J2, J3);
        err   = norm(Y(:) - Y_hat(:)) / norm_Y;
        fit_history(m) = err;
        iter_count = m;

        % Convergence check
        if abs(prev_err - err) / max(1e-10, prev_err) < eps_tol
            break;
        end
        if err < 1e-8
            break;
        end
        prev_err = err;
    end

    fit_history = fit_history(1:iter_count);
end


function KR = kr(A, B)
% Khatri-Rao product: KR(:,l) = kron(a_l, b_l)
% A: m×L,  B: n×L  →  KR: (m*n)×L
    [m, L] = size(A);
    n = size(B, 1);
    KR = zeros(m*n, L);
    for l = 1:L
        KR(:,l) = kron(A(:,l), B(:,l));
    end
end


function U = svd_init(M, L)
    % SVD-based initialization for better convergence
    % Use top L left singular vectors from SVD of mode unfolding
    if size(M, 1) > L
        [U, ~, ~] = svds(M, L);  % Top L singular vectors
    else
        % If dimension < L, use random init (shouldn't happen for tensors)
        [Q, ~] = qr(randn(size(M, 1), L) + 1j*randn(size(M, 1), L), 0);
        U = Q(:, 1:min(L, size(Q, 2)));
    end
    
    % Ensure normalization
    for l = 1:size(U, 2)
        U(:,l) = U(:,l) / (norm(U(:,l)) + eps);
    end
end


function Y_hat = reconstruct(A, B, C, J1, J2, J3)
% Y_(1) = A * (C⊙B)^T => reshape back to tensor
    CB = kr(C, B);  % (J3*J2) x L
    Y1_hat = A * CB.';  % J1 x (J3*J2)
    % Y1 was reshaped from (J1, J2, J3) where J2 varies fastest, then J3
    % Since reshape(Y, J1, J2*J3): columns go J2 fast, J3 slow
    % And kr(C,B) produces kron(c_l, b_l) = J3*J2 with b_l (J2) varying fast
    % This matches: Y1_hat has correct ordering
    Y_hat = reshape(Y1_hat, J1, J2, J3);
end
