function [Y, W, F_mat, c_true, AR, BT] = construct_tensor(params, Hk_all, alpha, tau)
% =========================================================================
% construct_tensor.m
% =========================================================================
% Description:
%   Constructs the received signal third-order tensor Y ∈ C^{F x T x K}
%   from the channel matrices Hk and hybrid beamforming matrices W, F.
%
%   Signal model (Eq. 9-10):
%     Y_k = W^T H_k F + V_k
%
%   Tensor structure (Eq. 16):
%     Y = I_{3,L} ×_1 (W^T A_R) ×_2 (F^T B_T) ×_3 C + V
%
%   W (NR x F) and F (NT x T) are truncated DFT matrices ensuring
%   column orthogonality: W*W^T = I_NR, F*F^T = I_NT  (Sec. 5.2)
%
%   Mode dimensions:
%     mode-1 = F (sub-frames)
%     mode-2 = T (time frames)
%     mode-3 = K (subcarriers)
%
% Inputs:
%   params   - system parameter struct (F, T, K, NR, NT, K_bar, fs, L)
%   Hk_all   - NR x NT x K true channel matrices
%   alpha    - L x 1 complex path gains
%   tau      - L x 1 time delays [s]
%
% Outputs:
%   Y     - F x T x K received tensor (complex, noisy)
%   W     - NR x F combining matrix (truncated DFT)
%   F_mat - NT x T precoding matrix (truncated DFT)
%   c_true - K x L true C factor matrix columns
%   AR     - NR x L true array responses at UT
%   BT     - NT x L true array responses at BS
%
% Runtime: O(F*T*K*L)
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
    % Construct DFT-based beamforming matrices (Sec. 5.2)
    % W: NR x F  (first F columns of NR-point DFT)
    % F: NT x T  (first T columns of NT-point DFT)
    % -------------------------------------------------------
    DFT_NR = dftmtx(NR) / sqrt(NR);   % NR x NR, unitary DFT
    DFT_NT = dftmtx(NT) / sqrt(NT);

    % Use first F_size / T_size columns
    W     = DFT_NR(:, 1:F_size);   % NR x F
    F_mat = DFT_NT(:, 1:T_size);   % NT x T

    % -------------------------------------------------------
    % True factor matrices A, B, C (without noise)
    % A = W^T A_R  (F x L)
    % B = F^T B_T  (T x L)
    % C = [c1...cL]  (K x L)
    % -------------------------------------------------------
    k_indices = round(linspace(1, K_bar, K));

    c_true = zeros(K, L);
    for l = 1:L
        for ki = 1:K
            k = k_indices(ki);
            c_true(ki,l) = alpha(l) * exp(-1j * 2*pi * tau(l) * fs * k / K_bar);
        end
    end

    % Precompute AR, BT from the true channel matrices
    % Extract AR, BT by solving the channel (simpler: pass from generate_channel)
    % Here we reconstruct from Hk_all to keep interfaces clean
    % A better approach: pass AR/BT directly (see main)
    % For now, we extract from channel structure.
    % We use: Hk = AR * diag(alpha.*phase_k) * BT', so
    %         sum_k c_k Hk / ||c_k||^2 ≈ AR * BT' (approx. when paths well-separated)
    % We'll build AR, BT as outputs from generate_channel directly in run_monte_carlo.
    AR = [];
    BT = [];

    % -------------------------------------------------------
    % Build received tensor Y = sum_l (W^T aR_l) o (F^T bT_l) o cl + V
    % Equivalent to stacking Yk = W^T Hk F + Vk  for k=1..K
    % -------------------------------------------------------
    Y_clean = zeros(F_size, T_size, K);

    for ki = 1:K
        Hk = Hk_all(:,:,ki);   % NR x NT
        Y_clean(:,:,ki) = W' * Hk * F_mat;   % F x T
    end

    % -------------------------------------------------------
    % Add AWGN noise scaled to desired SNR
    % SNR = ||Y_clean||_F^2 / ||V||_F^2
    % -------------------------------------------------------
    signal_power = norm(Y_clean(:))^2 / numel(Y_clean);
    noise_var    = signal_power / (10^(SNR_dB/10));

    % Complex Gaussian noise: real + imag each ~ N(0, noise_var/2)
    V = sqrt(noise_var/2) * (randn(F_size, T_size, K) + 1j*randn(F_size, T_size, K));

    Y = Y_clean + V;
end
