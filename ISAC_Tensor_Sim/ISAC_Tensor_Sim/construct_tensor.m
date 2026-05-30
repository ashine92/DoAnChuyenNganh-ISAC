function [Y, W, F_mat, c_true, AR, BT] = construct_tensor(params, Hk_all, alpha, tau, pl_all)
% =========================================================================
% construct_tensor.m
% =========================================================================
% Description:
%   Constructs the received signal tensor Y ∈ C^{F x T x K} and returns
%   the true beamforming matrices W, F and the true factor matrices.
%
%   Signal model (Eq. 10):  Y_k = W^T H_k F + V_k
%   Tensor CP form (Eq. 16): Y = Σ_l (W^T a_{R,l}) ∘ (F^T b_{T,l}) ∘ c_l + V
%
%   W (NR×F) and F_mat (NT×T): first F/T columns of the normalized DFT matrix.
%   Column orthogonality: W'*W = I_F,  F_mat'*F_mat = I_T.
%
% Inputs:
%   params   - struct: F, T, K, K_bar, NR, NT, fs, L, SNR_dB
%   Hk_all   - NR x NT x K channel matrices
%   alpha    - L x 1 complex gains
%   tau      - L x 1 time delays [s]
%   pl_all   - 3 x L SP positions [m]  (used to build AR, BT)
%
% Outputs:
%   Y      - F x T x K noisy received tensor
%   W      - NR x F combining matrix (truncated DFT)
%   F_mat  - NT x T precoding matrix (truncated DFT)
%   c_true - K x L true C-factor columns  c_l = alpha_l * c_bar(tau_l)
%   AR     - NR x L true array responses at UT
%   BT     - NT x L true array responses at BS
%
% Runtime: O(F*T*K + NR*NT*K)
% =========================================================================

    F_size = params.F;
    T_size = params.T;
    K      = params.K;
    K_bar  = params.K_bar;
    NR     = params.NR;
    NT     = params.NT;
    fs     = params.fs;
    L      = params.L;
    SNR_dB = params.SNR_dB;

    % -------------------------------------------------------
    % DFT beamforming matrices (Sec 5.2, paper)
    % W'*W = I_F  and  F_mat'*F_mat = I_T  by unitary DFT property
    % (Note: If F > NR or T > NT, we use an oversampled DFT dictionary)
    % -------------------------------------------------------
    if F_size <= NR
        W = dftmtx(NR) / sqrt(NR);
        W = W(:, 1:F_size);   % NR x F
    else
        W_full = dftmtx(F_size) / sqrt(F_size);
        W = W_full(1:NR, :);  % NR x F
    end

    if T_size <= NT
        F_mat = dftmtx(NT) / sqrt(NT);
        F_mat = F_mat(:, 1:T_size); % NT x T
    else
        F_full = dftmtx(T_size) / sqrt(T_size);
        F_mat = F_full(1:NT, :);    % NT x T
    end

    % -------------------------------------------------------
    % Subcarrier index set (K uniformly selected from K_bar)
    % -------------------------------------------------------
    k_indices = round(linspace(1, K_bar, K));   % 1 x K

    % -------------------------------------------------------
    % True C-factor matrix: c_true(:,l) = alpha_l * c_bar(tau_l)
    %   c_bar(tau, k) = exp(-j2*pi*tau*fs*k/K_bar)
    % -------------------------------------------------------
    c_true = zeros(K, L);
    for l = 1:L
        c_true(:,l) = alpha(l) * exp(-1j * 2*pi * tau(l) * fs * k_indices' / K_bar);
    end

    % -------------------------------------------------------
    % True AR, BT from near-field array geometry
    % -------------------------------------------------------
    AR = zeros(NR, L);
    BT = zeros(NT, L);
    if nargin >= 5 && ~isempty(pl_all)
        for l = 1:L
            [AR(:,l), BT(:,l)] = near_field_array_response(params, pl_all(:,l));
        end
    end

    % -------------------------------------------------------
    % Build noiseless received tensor: Y_clean(:,:,k) = W' * Hk * F_mat
    % -------------------------------------------------------
    Y_clean = zeros(F_size, T_size, K);
    for ki = 1:K
        Y_clean(:,:,ki) = W' * Hk_all(:,:,ki) * F_mat;   % F x T
    end

    % -------------------------------------------------------
    % Add AWGN: SNR = ||Y_clean||_F^2 / ||V||_F^2
    % -------------------------------------------------------
    sig_pwr   = norm(Y_clean(:))^2 / numel(Y_clean);
    noise_var = sig_pwr / 10^(SNR_dB/10);
    V = sqrt(noise_var/2) * (randn(F_size,T_size,K) + 1j*randn(F_size,T_size,K));

    Y = Y_clean + V;
end
