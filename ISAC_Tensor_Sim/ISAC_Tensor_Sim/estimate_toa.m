function tau_hat = estimate_toa(C_hat, params)
% =========================================================================
% estimate_toa.m  - Robust 1D correlation search
% =========================================================================
% * Mathematical Background:
%   Estimates Time of Arrival (ToA) from C_hat by maximizing correlation:
%   tau_hat = argmax_tau | c_hat^H * c_true(tau) |
%
%   This avoids the phase unwrapping ambiguity that occurs when K is small
%   and adjacent subcarrier phase steps exceed pi.
% =========================================================================

    K     = params.K;
    K_bar = params.K_bar;
    fs    = params.fs;
    L     = params.L;

    k_indices = round(linspace(1, K_bar, K))';
    tau_hat = zeros(L, 1);

    % Search grid: 0 to 200 ns
    % For higher precision, we do a two-step search (coarse then fine)
    tau_grid_coarse = linspace(0, 200e-9, 1000); 

    for l = 1:L
        cl_hat = C_hat(:, l);
        
        % Coarse search
        C_dict_coarse = exp(-1j * 2*pi * fs * k_indices * tau_grid_coarse / K_bar); % K x 1000
        corr_coarse = abs(cl_hat' * C_dict_coarse);
        [~, best_idx_c] = max(corr_coarse);
        tau_c = tau_grid_coarse(best_idx_c);
        
        % Fine search around the coarse estimate (+/- 1 coarse step)
        step_c = tau_grid_coarse(2) - tau_grid_coarse(1);
        tau_grid_fine = linspace(max(0, tau_c - step_c), tau_c + step_c, 1000);
        
        C_dict_fine = exp(-1j * 2*pi * fs * k_indices * tau_grid_fine / K_bar);
        corr_fine = abs(cl_hat' * C_dict_fine).^2;
        [~, best_idx_f] = max(corr_fine);
        
        % Off-grid Parabolic Interpolation to break the grid resolution limit
        if best_idx_f > 1 && best_idx_f < length(tau_grid_fine)
            y1 = corr_fine(best_idx_f - 1);
            y2 = corr_fine(best_idx_f);
            y3 = corr_fine(best_idx_f + 1);
            delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3 + eps);
            step_f = tau_grid_fine(2) - tau_grid_fine(1);
            tau_hat(l) = tau_grid_fine(best_idx_f) + delta * step_f;
        else
            tau_hat(l) = tau_grid_fine(best_idx_f);
        end
    end
end
