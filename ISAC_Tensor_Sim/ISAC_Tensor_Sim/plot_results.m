% =========================================================================
% plot_results.m
% =========================================================================
% Description:
%   Collection of plot functions for reproducing Figures 4-8 from the paper.
%   All figures are saved as PNG files in the current directory.
%
%   Functions:
%     plot_fig4(results, params)   - NMSE vs SNR (6 subplots)
%     plot_fig5(results, params)   - NMSE vs K   (6 subplots)
%     plot_fig6(results, cases, p) - NMSE angles vs SNR (4 subplots)
%     plot_fig7(results, cases, p) - NMSE ToA/Loc vs SNR (3 subplots)
%     plot_fig8(results, cases, p) - 3D localization (4 subplots)
%
% =========================================================================

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


function plot_fig8(results_cell, cases, params)
% 3D localization visualization (4 subplots)

    case_labels = {'(a) L=3, F=T=50','(b) L=3, F=T=80', ...
                   '(c) L=2, F=T=50','(d) L=2, F=T=80'};

    fig = figure('Position',[50 50 1000 900], 'Name', 'Figure 8');

    % Shared legend items
    leg_entries = {};

    for ci = 1:4
        r = results_cell{ci};
        ax = subplot(2,2,ci);
        hold on; grid on; box on;
        view(3);

        % BS (known, red square)
        pT = r.pT;
        h1 = plot3(pT(1), pT(2), pT(3), 'rs', 'MarkerSize', 10, ...
            'MarkerFaceColor','r', 'DisplayName','BS');

        % UT true (circle)
        pR = r.pR_true;
        h2 = plot3(pR(1), pR(2), pR(3), 'ko', 'MarkerSize', 8, ...
            'MarkerFaceColor','k', 'DisplayName','UT true');

        % UT estimated (circle with 'Proposed')
        pR_h = r.pR_hat;
        h3 = plot3(pR_h(1), pR_h(2), pR_h(3), 'g^', 'MarkerSize', 8, ...
            'MarkerFaceColor','g', 'DisplayName','UT est');

        % Connect true and estimated UT
        plot3([pR(1),pR_h(1)],[pR(2),pR_h(2)],[pR(3),pR_h(3)], 'g--','LineWidth',1);

        % SPs true (blue squares)
        L_ci = cases(ci).L;
        for l = 1:L_ci
            pl_t = r.pl_true(:,l);
            h4 = plot3(pl_t(1), pl_t(2), pl_t(3), 'bs', 'MarkerSize', 8, ...
                'MarkerFaceColor','b', 'DisplayName','SPs true');
        end

        % SPs estimated
        for l = 1:L_ci
            pl_e = r.pl_hat(:,l);
            h5 = plot3(pl_e(1), pl_e(2), pl_e(3), 'ms', 'MarkerSize', 8, ...
                'MarkerFaceColor','m', 'DisplayName','SPs est');

            % Connect true and estimated SP
            plot3([r.pl_true(1,l),pl_e(1)],[r.pl_true(2,l),pl_e(2)], ...
                  [r.pl_true(3,l),pl_e(3)], 'm--','LineWidth',1);

            % Draw lines from BS to SP and SP to UT (signal path)
            plot3([pT(1),r.pl_true(1,l)],[pT(2),r.pl_true(2,l)],[pT(3),r.pl_true(3,l)],'b-','LineWidth',0.5,'HandleVisibility','off');
            plot3([r.pl_true(1,l),pR(1)],[r.pl_true(2,l),pR(2)],[r.pl_true(3,l),pR(3)],'b-','LineWidth',0.5,'HandleVisibility','off');
        end

        xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
        title(sprintf('%s', case_labels{ci}), 'FontSize', 9);

        % Only add legend to first subplot to avoid repetition
        if ci == 1
            legend([h1,h2,h3,h4,h5], {'BS','UT true','UT est','SPs true','SPs est'}, ...
                'Location','best','FontSize',7);
        end

        set(gca,'FontSize',8);
    end

    sgtitle(sprintf('Fig.8: 3D Localization Visualization (SNR=%d dB, K=%d)', ...
        params.SNR_dB, params.K), 'FontSize', 11);
    saveas(fig, 'Figure8_3D_Localization.png');
    fprintf('  Saved: Figure8_3D_Localization.png\n');
end
