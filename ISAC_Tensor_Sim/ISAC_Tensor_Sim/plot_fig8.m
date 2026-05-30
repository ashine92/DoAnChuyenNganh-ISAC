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
