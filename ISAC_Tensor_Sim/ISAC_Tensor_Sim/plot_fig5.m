function plot_fig5(results, params)
% NMSE vs K for different algorithms

    K_vec   = results.K_vec;
    fields  = {'az_R','el_R','az_T','el_T','pl','pR'};
    ylabels = {'NMSE (\theta^{az}_R)', 'NMSE (\theta^{el}_R)', ...
               'NMSE (\theta^{az}_T)', 'NMSE (\theta^{el}_T)', ...
               'NMSE (p_l)', 'NMSE (p_R)'};
    panel_labels = {'(a) AoA azimuth','(b) AoA elevation', ...
                    '(c) AoD azimuth','(d) AoD elevation', ...
                    '(e) SPs position','(f) UT position'};

    colors  = {[0.8 0 0], [0 0.5 0], [0 0 0.8], [0.5 0 0.5]};
    markers = {'o-','s--','^-','d-.'};
    methods = {'Proposed','MUSIC_LSPS','PUDD','CRB'};
    leg_str = {'Proposed','MUSIC-LSPS','PUDD','CRB'};

    fig = figure('Position',[100 100 1100 700], 'Name', 'Figure 5');
    for pi_ = 1:6
        subplot(2,3,pi_);
        hold on; grid on; box on;
        for mi = 1:3
            meth_data = results.(methods{mi}).(fields{pi_});
            semilogy(K_vec, meth_data, markers{mi}, 'Color', colors{mi}, ...
                'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', leg_str{mi});
        end
        crb_data = results.CRB.(fields{pi_});
        semilogy(K_vec, crb_data, markers{4}, 'Color', colors{4}, ...
            'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', leg_str{4});
        xlabel('K', 'FontSize', 10);
        ylabel(ylabels{pi_}, 'FontSize', 10);
        title(panel_labels{pi_}, 'FontSize', 9);
        legend('Location','northeast', 'FontSize', 8);
        xlim([K_vec(1), K_vec(end)]);
        set(gca, 'FontSize', 9);
    end
    sgtitle(sprintf('Fig.5: NMSE vs K  (F=T=%d, SNR=%d dB, L=%d)', ...
        params.F, params.SNR_dB, params.L), 'FontSize', 11);
    saveas(fig, 'Figure5_NMSE_vs_K.png');
    fprintf('  Saved: Figure5_NMSE_vs_K.png\n');
end
