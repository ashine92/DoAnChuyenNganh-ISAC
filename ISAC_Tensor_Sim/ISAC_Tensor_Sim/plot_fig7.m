function plot_fig7(results_cell, cases, params)
% NMSE of ToA and localization vs SNR for 4 cases

    SNR_vec = results_cell{1}.SNR_vec;
    fields  = {'tau','pR','pl'};
    ylabels = {'NMSE (\tau)', 'NMSE (p_R)', 'NMSE (p_l)'};

    colors  = {[0.8 0 0], [0.8 0 0], [0 0 0.8], [0 0 0.8]};
    styles  = {'-o','-s','--o','--s'};

    fig = figure('Position',[100 100 900 350], 'Name', 'Figure 7');
    for pi_ = 1:3
        subplot(1,3,pi_);
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
    sgtitle(sprintf('Fig.7: ToA & Localization NMSE vs SNR  (K=%d)', params.K), 'FontSize', 11);
    saveas(fig, 'Figure7_Localization_NMSE_vs_SNR.png');
    fprintf('  Saved: Figure7_Localization_NMSE_vs_SNR.png\n');
end
