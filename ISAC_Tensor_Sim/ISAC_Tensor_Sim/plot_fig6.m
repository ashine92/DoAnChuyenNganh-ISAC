function plot_fig6(results, params)
% =========================================================================
% plot_fig6.m — Plot RMSE vs N (Figure 6)
% =========================================================================

    N_vec   = results.N_vec.^2;
    fields  = {'az_R','el_R','az_T','el_T','pl','pR'};
    ylabels = {'RMSE (\theta^{az}_R)', 'RMSE (\theta^{el}_R)', ...
               'RMSE (\theta^{az}_T)', 'RMSE (\theta^{el}_T)', ...
               'RMSE (p_l)', 'RMSE (p_R)'};
    panel_labels = {'(a) AoA azimuth','(b) AoA elevation', ...
                    '(c) AoD azimuth','(d) AoD elevation', ...
                    '(e) SPs position','(f) UT position'};

    colors  = {[0.8 0 0], [0 0.5 0], [0 0 0.8], [0.5 0 0.5]};
    markers = {'o-','s--','^-','d-.'};
    methods = {'Proposed','MUSIC_LSPS','PUDD','CRB'};
    leg_str = {'Proposed','MUSIC-LSPS','PUDD','CRB'};

    fig = figure('Position',[100 100 1100 700], 'Name', 'Figure 6');
    for pi_ = 1:6
        subplot(2,3,pi_);
        hold on; grid on; box on;
        for mi = 1:3
            if isfield(results, methods{mi})
                data = results.(methods{mi}).(fields{pi_});
                if any(data > 0)
                    data = max(data, 1e-12);
                    semilogy(N_vec, data, markers{mi}, 'Color', colors{mi}, ...
                        'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', leg_str{mi});
                end
            end
        end
        if isfield(results, 'CRB')
            crb_data = results.CRB.(fields{pi_});
            if any(crb_data > 0)
                crb_data = max(crb_data, 1e-12);
                semilogy(N_vec, crb_data, markers{4}, 'Color', colors{4}, ...
                    'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', leg_str{4});
            end
        end
        xlabel('N', 'FontSize', 10);
        ylabel(ylabels{pi_}, 'FontSize', 10);
        title(panel_labels{pi_}, 'FontSize', 9);
        legend('Location','southwest', 'FontSize', 7);
        xlim([N_vec(1), N_vec(end)]);
        set(gca, 'FontSize', 9, 'YScale', 'log');
    end
    sgtitle(sprintf('Fig.6: RMSE vs N  (K=%d, SNR=%d dB, L=%d)', ...
        params.K, params.SNR_dB, params.L), 'FontSize', 11);
    saveas(fig, 'Figure6_RMSE_vs_N.png');
    fprintf('  Saved: Figure6_RMSE_vs_N.png\n');
end
