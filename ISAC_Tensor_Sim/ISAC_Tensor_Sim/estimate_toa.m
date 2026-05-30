function tau_hat = estimate_toa(C_hat, params)
% =========================================================================
% estimate_toa.m
% =========================================================================
% Description:
%   Estimates the Time of Arrival (ToA) for each path from the estimated
%   factor matrix C_hat using maximum likelihood (column correlation).
%
%   Mathematical derivation (Appendix A, Eq. 32):
%     tau_hat_l = argmax_{tau_l} |c_hat_l^H * c_bar(tau_l)|^2
%                                 / (||c_bar(tau_l)||^2 * ||c_hat_l||^2)
%
%   where c_bar(tau) = [exp(-j2*pi*tau*fs*k1/K_bar), ...,
%                        exp(-j2*pi*tau*fs*kK/K_bar)]^T  is the steering
%   vector for a candidate delay tau.
%
%   One-dimensional exhaustive search over Ns uniformly spaced points
%   in the range [0, K_bar/fs].
%
% Inputs:
%   C_hat  - K x L estimated factor matrix (subcarrier dimension)
%   params - struct with fields: K, K_bar, fs, Ns, L
%
% Outputs:
%   tau_hat - L x 1 estimated time delays [s]
%
% Runtime: O(Ns * L * K)
% =========================================================================

    K     = params.K;
    K_bar = params.K_bar;
    fs    = params.fs;
    Ns    = params.Ns;
    L     = params.L;

    % Search range: tau in [0, K_bar/fs]
    tau_max    = K_bar / fs;
    tau_search = linspace(0, tau_max, Ns);

    % Subcarrier indices (same selection as in channel generation)
    k_indices = round(linspace(1, K_bar, K))';   % K x 1

    % Precompute steering matrix: K x Ns
    % c_bar(tau, k) = exp(-j2*pi*tau*fs*k/K_bar)
    % Phase matrix: K x Ns
    phase_mat = exp(-1j * 2*pi/K_bar * fs * (k_indices * tau_search));  % K x Ns

    tau_hat = zeros(L, 1);

    for l = 1:L
        cl_hat = C_hat(:, l);              % K x 1
        cl_hat_norm = norm(cl_hat);

        % Correlation with each candidate: 1 x Ns
        corr_vals = abs(cl_hat' * phase_mat).^2;

        % Normalize by ||c_bar||^2 * ||c_hat||^2
        c_bar_norm_sq = sum(abs(phase_mat).^2, 1);   % 1 x Ns (all equal to K)
        corr_norm     = corr_vals ./ (c_bar_norm_sq * cl_hat_norm^2 + eps);

        [~, idx]   = max(corr_norm);
        tau_hat(l) = tau_search(idx);
    end
end
