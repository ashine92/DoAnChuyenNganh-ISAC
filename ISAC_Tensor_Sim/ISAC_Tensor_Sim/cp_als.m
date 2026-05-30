function [A_hat, B_hat, C_hat, iter_count, fit_history] = cp_als(Y, L, params)
% =========================================================================
% cp_als.m
% =========================================================================
% Description:
%   Computes the CP (PARAFAC/CANDECOMP) decomposition of a third-order
%   tensor Y ∈ C^{J1 x J2 x J3} into L rank-1 components using the
%   Alternating Least Squares (ALS) algorithm.
%
%   Mathematical model (Eq. 28):
%     Y ≈ sum_{l=1}^{L} a_l ∘ b_l ∘ c_l
%
%   Mode unfoldings used (Eq. 25):
%     [Y]_(1)^T = A * (C ⊙ B)^T
%     [Y]_(2)^T = B * (C ⊙ A)^T
%     [Y]_(3)^T = C * (B ⊙ A)^T
%
%   ALS updates (Eq. 29-30):
%     A ← [Y]_(1) * (C⊙B) * [(C^T C * B^T B)]^{-1}
%     B ← [Y]_(2) * (C⊙A) * [(C^T C * A^T A)]^{-1}
%     C ← [Y]_(3) * (B⊙A) * [(B^T B * A^T A)]^{-1}
%
%   Initialization: SVD of each mode unfolding (Eq. 29 initialization).
%
% Inputs:
%   Y      - J1 x J2 x J3 third-order complex tensor
%   L      - number of CP components (rank)
%   params - struct with fields: max_iter, eps_tol
%
% Outputs:
%   A_hat       - J1 x L factor matrix (sub-frames dimension)
%   B_hat       - J2 x L factor matrix (time frames dimension)
%   C_hat       - J3 x L factor matrix (subcarrier dimension)
%   iter_count  - number of ALS iterations performed
%   fit_history - vector of ||Y - Y_hat||_F / ||Y||_F at each iteration
%
% Khatri-Rao product: khatri_rao(A,B) defined in utils/khatri_rao.m
%
% Runtime: O(max_iter * (J1*J2*J3*L + J1*J3*L^2 + L^3))
%          ≈ O(FTKL) per iteration
% =========================================================================

    max_iter = params.max_iter;
    eps_tol  = params.eps_tol;

    [J1, J2, J3] = size(Y);

    % -------------------------------------------------------
    % Compute mode-n unfoldings of Y
    % Convention: [Y]_(n) has mode-n fibers as columns
    % For a (J1 x J2 x J3) tensor:
    %   [Y]_(1): J1 x (J2*J3)  — rows are mode-1 fibers
    %   [Y]_(2): J2 x (J1*J3)
    %   [Y]_(3): J3 x (J1*J2)
    %
    % Paper uses transposed convention (Eq. 5): [Y]_(n)^T
    % We compute [Y]_(n) as the standard unfolding.
    % -------------------------------------------------------
    Y1 = tensor_unfold(Y, 1);   % J1 x (J2*J3) — standard mode-1
    Y2 = tensor_unfold(Y, 2);   % J2 x (J1*J3)
    Y3 = tensor_unfold(Y, 3);   % J3 x (J1*J2)

    norm_Y = norm(Y(:));

    % -------------------------------------------------------
    % Initialization: first L left singular vectors of each unfolding
    % -------------------------------------------------------
    A_hat = svd_init(Y1, L);
    B_hat = svd_init(Y2, L);
    C_hat = svd_init(Y3, L);

    fit_history  = zeros(max_iter, 1);
    iter_count   = 0;
    prev_fit     = inf;

    % -------------------------------------------------------
    % ALS main loop
    % -------------------------------------------------------
    for m = 1:max_iter
        % --- Update A ---
        % A ← Y_(1) * (C⊙B) * [(C^T C * B^T B)]^{-1}
        % Note: Y_(1) is J1 x (J2*J3);  C⊙B is (J2*J3) x L
        CB = khatri_rao(C_hat, B_hat);          % (J3*J2) x L  [or J2*J3]
        % Need ordering: [Y]_(1) = A*(C⊙B)^T  means columns of (C⊙B) index mode2,3
        % Standard: Y1 = J1 x (J2*J3), unfold order: j3 fastest then j2
        % We use consistent ordering throughout
        gram = (C_hat' * C_hat) .* (B_hat' * B_hat);   % L x L Hadamard
        A_hat = Y1 * CB / gram;   % J1 x L  (least squares)

        % --- Update B ---
        CA = khatri_rao(C_hat, A_hat);           % (J3*J1) x L
        gram = (C_hat' * C_hat) .* (A_hat' * A_hat);
        B_hat = Y2 * CA / gram;   % J2 x L

        % --- Update C ---
        BA = khatri_rao(B_hat, A_hat);           % (J2*J1) x L
        gram = (B_hat' * B_hat) .* (A_hat' * A_hat);
        C_hat = Y3 * BA / gram;   % J3 x L

        % --- Compute fit ---
        Y_hat = cp_reconstruct(A_hat, B_hat, C_hat);
        residual = norm(Y(:) - Y_hat(:)) / norm_Y;
        fit_history(m) = residual;
        iter_count = m;

        % --- Convergence check ---
        if abs(prev_fit - residual) < eps_tol && m > 5
            break;
        end
        prev_fit = residual;
    end

    fit_history = fit_history(1:iter_count);
end


%% --- Helper: Tensor unfolding ---
function U = tensor_unfold(X, mode)
% Mode-n unfolding of tensor X.
% For (J1 x J2 x J3):
%   mode=1: U is J1 x (J2*J3), columns are mode-1 fibers
%           arranged so that index j2 varies faster than j3
%   mode=2: U is J2 x (J1*J3), j1 varies faster than j3
%   mode=3: U is J3 x (J1*J2), j1 varies faster than j2
%
% This ordering matches the Khatri-Rao product convention used in ALS.

    sz = size(X);
    ndims_X = length(sz);
    % Permute so that `mode` is the first index, then reshape
    order = [mode, setdiff(1:ndims_X, mode)];
    X_perm = permute(X, order);
    U = reshape(X_perm, sz(mode), []);
end


%% --- Helper: SVD initialization ---
function U = svd_init(M, L)
% Returns first L left singular vectors of matrix M.
    [U, ~, ~] = svds(M, L);
    if size(U,2) < L
        % Pad with random columns if svds returns fewer
        pad = randn(size(U,1), L-size(U,2)) + 1j*randn(size(U,1), L-size(U,2));
        U = [U, pad];
    end
end


%% --- Helper: CP reconstruction ---
function Y_hat = cp_reconstruct(A, B, C)
% Reconstructs tensor from factor matrices A (J1xL), B (J2xL), C (J3xL).
% Y_hat(j1,j2,j3) = sum_l A(j1,l)*B(j2,l)*C(j3,l)

    [J1, L] = size(A);
    J2 = size(B, 1);
    J3 = size(C, 1);

    % Vectorized outer product construction
    % Y_(3) = C * (B⊙A)^T  → reshape to tensor
    BA = khatri_rao(B, A);    % (J2*J1) x L
    Y3_hat = C * BA';         % J3 x (J2*J1)

    % Reshape: Y3_hat is [J3 x (J2*J1)] with j1 fast, j2 slow
    % We need [J1 x J2 x J3]
    % Y_(3) has ordering j1 fastest — matches tensor_unfold(mode=3)
    Y_hat = permute(reshape(Y3_hat', [J1, J2, J3]), [1, 2, 3]);
end
