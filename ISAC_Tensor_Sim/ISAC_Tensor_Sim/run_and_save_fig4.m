% Run Figure 4 and save key metrics to file for analysis
clear all; close all; clc;

% Run main_fig4
try
    main_fig4;
    
    % The results should be in the workspace
    % Try to extract and save metrics
    
    fid = fopen('fig4_results_latest.txt', 'w');
    fprintf(fid, 'Figure 4 Execution Results\n');
    fprintf(fid, '=========================\n\n');
    
    % Check if NMSE_proposed exists
    if exist('NMSE_proposed', 'var')
        fprintf(fid, 'NMSE_proposed shape: %d x %d\n', size(NMSE_proposed, 1), size(NMSE_proposed, 2));
        fprintf(fid, 'NMSE_proposed (first row): ');
        fprintf(fid, '%e ', NMSE_proposed(1, :));
        fprintf(fid, '\n\n');
    end
    
    if exist('NMSE_music', 'var')
        fprintf(fid, 'NMSE_music (first row): ');
        fprintf(fid, '%e ', NMSE_music(1, :));
        fprintf(fid, '\n\n');
    end
    
    if exist('SNR_dB', 'var')
        fprintf(fid, 'SNR_dB: ');
        fprintf(fid, '%d ', SNR_dB);
        fprintf(fid, '\n\n');
    end
    
    % Check localization results
    if exist('pR_nmse', 'var')
        fprintf(fid, 'pR_nmse (UT localization): ');
        fprintf(fid, '%e ', pR_nmse);
        fprintf(fid, '\n\n');
    end
    
    fprintf(fid, 'Execution completed successfully\n');
    fclose(fid);
    
    disp('Results saved to fig4_results_latest.txt');
    
catch ME
    fid = fopen('fig4_results_latest.txt', 'w');
    fprintf(fid, 'ERROR: %s\n', ME.message);
    fprintf(fid, 'Stack:\n');
    fprintf(fid, '%s\n', ME.stack);
    fclose(fid);
    disp(ME.message);
end

exit;
