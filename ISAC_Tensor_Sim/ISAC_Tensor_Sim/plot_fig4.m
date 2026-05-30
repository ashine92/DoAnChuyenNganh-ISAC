function plot_fig4(results, params)
% NMSE vs SNR for proposed, MUSIC-LSPS, PUDD, CRB
% Reproduces Figure 4: 6 subplots (az_R, el_R, az_T, el_T, SPs pos, UT pos)

    SNR_vec = results.SNR_vec;
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

    fig = figure('Position',[100 100 1100 700], 'Name', 'Figure 4');
    for pi_ = 1:6
        subplot(2,3,pi_);
        hold on; grid on; box on;

        for mi = 1:3
            meth_data = results.(methods{mi}).(fields{pi_});
            semilogy(SNR_vec, meth_data, markers{mi}, 'Color', colors{mi}, ...
                'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', leg_str{mi});
        end
        % CRB
        crb_data = results.CRB.(fields{pi_});
        semilogy(SNR_vec, crb_data, markers{4}, 'Color', colors{4}, ...
            'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', leg_str{4});

        xlabel('SNR (dB)', 'FontSize', 10);
        ylabel(ylabels{pi_}, 'FontSize', 10);
        title(panel_labels{pi_}, 'FontSize', 9);
        legend('Location','northeast', 'FontSize', 8);
        xlim([SNR_vec(1), SNR_vec(end)]);
        set(gca, 'FontSize', 9);
    end

    sgtitle(sprintf('Fig.4: NMSE vs SNR  (F=T=%d, K=%d, L=%d)', ...
        params.F, params.K, params.L), 'FontSize', 11);

    saveas(fig, 'Figure4_NMSE_vs_SNR.png');
    fprintf('  Saved: Figure4_NMSE_vs_SNR.png\n');
end
