function plot_fig6(results_cell, cases, params)
% NMSE of angle parameters vs SNR for 4 cases

    SNR_vec = results_cell{1}.SNR_vec;
    fields  = {'az_R','el_R','az_T','el_T'};
    ylabels = {'NMSE (\theta^{az}_R)', 'NMSE (\theta^{el}_R)', ...
               'NMSE (\theta^{az}_T)', 'NMSE (\theta^{el}_T)'};

    colors  = {[0.8 0 0], [0.8 0 0], [0 0 0.8], [0 0 0.8]};
    styles  = {'-o','-s','--o','--s'};

    fig = figure('Position',[100 100 1000 500], 'Name', 'Figure 6');
    for pi_ = 1:4
        subplot(1,4,pi_);
        hold on; grid on; box on;
        for ci = 1:4
            data = results_cell{ci}.Proposed.(fields{pi_});
            semilogy(SNR_vec, data, styles{ci}, 'Color', colors{ci}, ...
                'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', cases(ci).label);
        end
        xlabel('SNR (dB)', 'FontSize', 10);
        ylabel(ylabels{pi_}, 'FontSize', 10);
        legend('Location','northeast', 'FontSize', 7);
        xlim([SNR_vec(1), SNR_vec(end)]);
        set(gca, 'FontSize', 9);
    end
    sgtitle(sprintf('Fig.6: Angle NMSE vs SNR  (K=%d)', params.K), 'FontSize', 11);
    saveas(fig, 'Figure6_Angle_NMSE_vs_SNR.png');
    fprintf('  Saved: Figure6_Angle_NMSE_vs_SNR.png\n');
end
