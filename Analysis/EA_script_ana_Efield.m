function EA_script_ana_Efield(EA_results)
%EA_SCRIPT_ANA_EFIELD - Plot electric field amplitude from ElectroAbsorption (EA)
% in a range of background light intensities or applied DC voltages
% this file is heavily based on the equivalent script for Impedance
% Spectroscopy: IS_script_ana_impedance
% The amplitude of the spatially averaged electric field first harmonic
% $(x_2-x_1)^-1 \int_{x_1}{x_2}EDC(x) * EAC(x) dx$ and second
% harmonic $(x_2-x_1)^-1 \int_{x_1}{x_2} EAC^2(x) dx$ are plotted
% versus the applied voltage frequency.
%
% Syntax:  EA_script_ana_Efield(EA_results)
%
% Inputs:
%   EA_RESULTS - a struct containing the most important results of the EA simulation
%
% Example:
%   EA_script_ana_Efield(EA_results)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also EA_script.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Davide Moia, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% March 2018; Last revision: March 2018

%------------- BEGIN CODE --------------

% check which was the variable being explored
if numel(unique(EA_results.Int)) > 1
    legend_text = EA_results.Int;
    legend_append = ' sun';
else
    legend_text = EA_results.Vdc;
    legend_append = ' Vdc';
end

% create a color array with one color more than necessary
jet_matrix = jet(length(legend_text) + 1);
% find the yellow (which in RGB code is 1,1,0) and remove it from the
% colors list
jet_yellow_logical = ismember(jet_matrix, [1, 1, 0], 'rows');
jet_no_yellow = jet_matrix(~jet_yellow_logical, :);
Int_colors = colormap(jet_no_yellow);

% round to two significant digits
legend_flip = round(legend_text, 2, 'significant');
% flip array and matrixes for starting from dark
legend_flip = string(flipud(legend_flip));
% add sun to numbers in legend
legend_flip = strcat(legend_flip, legend_append);
% replace zero in legend with dark
legend_flip(legend_flip=="0 sun") = "dark";

% preallocate figures handles
h = zeros(length(legend_text), 1);

%% plot
figure('Name', 'Amplitude of EA first harmonic E AC x E DC', 'NumberTitle', 'off');
    hold off
    for i = 1:length(legend_text)
        h(i) = plot(EA_results.Freq(i, :), EA_results.AC_ExDC_E_amp(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        if ~isnan(EA_results.AC_ExDC_E_i_amp(1))
            % due to ions
            plot(EA_results.Freq(i, :), EA_results.AC_ExDC_E_i_amp(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--');
        end
    end
    ax = gca;
    ax.XScale = 'log'; % for putting the scale in log
    ax.YScale = 'log'; % for putting the scale in log
    xlim([min(min(EA_results.Freq)), max(max(EA_results.Freq))])
    xlabel('Frequency [Hz]');
    ylabel('Abs(E_{AC} x E_{DC}) [V^2/cm^2]');
    legend(flipud(h), legend_flip)
    legend boxoff

figure('Name', 'Amplitude of EA second harmonic E_{AC}^2', 'NumberTitle', 'off');
    hold off
    for i = 1:length(legend_text)
        h(i) = plot(EA_results.Freq(i, :), EA_results.AC_Efield2_amp(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        plot(EA_results.Freq(i, :), EA_results.AC_Efield_amp_squared_mean(i, :)',...
            'Color', Int_colors(i, :), 'LineStyle', '-.');
        if ~isnan(EA_results.AC_Efield2_i_amp(1))
            % due to ions
            plot(EA_results.Freq(i, :), EA_results.AC_Efield2_i_amp(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--');
        end
    end
    ax = gca;
    ax.XScale = 'log'; % for putting the scale in log
    ax.YScale = 'log'; % for putting the scale in log
    xlim([min(min(EA_results.Freq)), max(max(EA_results.Freq))])
    xlabel('Frequency [Hz]');
    ylabel('Abs(E_{AC}^2) [V^2/cm^2]');
    legend(flipud(h), legend_flip)
    legend boxoff
end

%------------- END OF CODE --------------
