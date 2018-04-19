function ISwave_full_analysis_phase(ISwave_struct)
%ISWAVE_FULL_ANALYSIS_PHASE - Plot phase Bode plots from Impedance Spectroscopy (IS)
% in a range of background light intensities or applied DC voltages
% This script is very similar to EA_full_analysis_phase
%
% Syntax:  ISwave_full_analysis_phase(ISwave_struct)
%
% Inputs:
%   ISWAVE_STRUCT - a struct containing the most important results of the ISwave simulation
%
% Example:
%   ISwave_full_analysis_phase(ISwave_struct)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also ISwave_full_exec, EA_full_analysis_phase, ISwave_full_analysis_nyquist, IS_full_analysis_impedance.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: March 2018

%------------- BEGIN CODE --------------

% increase graphics font size
set(0, 'defaultAxesFontSize', 24);
% set image dimension
set(0, 'defaultfigureposition', [0, 0, 1000, 750]);
% set line thickness
set(0, 'defaultLineLineWidth', 2);

% check which was the variable being explored
if numel(unique(ISwave_struct.Int)) > 1
    legend_text = ISwave_struct.Int;
    legend_append = ' sun';
else
    legend_text = ISwave_struct.Vdc;
    legend_append = ' Vdc';
end

% create a color array with one color more than necessary
jet_matrix = jet(length(legend_text) + 1);
% find the yellow (which in RGB code is 1,1,0) and remove it from the
% colors list
jet_yellow_logical = ismember(jet_matrix, [1, 1, 0], 'rows');
jet_no_yellow = jet_matrix(~jet_yellow_logical, :);
jet_no_yellow_flip = flipud(jet_no_yellow);
Int_colors = colormap(jet_no_yellow_flip);

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

% in case just a single frequency was simulated, min and max are
% coincident, so they are modified by a small percentage
xlim_array = [min(min(ISwave_struct.Freq))*0.99, max(max(ISwave_struct.Freq)*1.01)];

%% do ISwave phase plot

% we want to plot phase of impedance, which is -phase(I vs V)
phase_n_deg = rad2deg(ISwave_struct.J_phase);
phase_i_deg = rad2deg(ISwave_struct.J_i_phase);
phase_U_deg = rad2deg(ISwave_struct.J_U_phase);
%phase_n_noionic_deg = rad2deg(results_struct.J_noionic_phase);
phase_dQ_deg = rad2deg(ISwave_struct.dQ_phase);

figure('Name', 'Phase plot of ISwave. Dashed: ionic; dotted: recombination; dashdotted: stored charge variation', 'NumberTitle', 'off')
    hold off
    for i = 1:length(legend_text)
        h(i) = plot(ISwave_struct.Freq(i, :), -phase_n_deg(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        plot(ISwave_struct.Freq(i, :), -phase_i_deg(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--', 'LineWidth', 1);
        plot(ISwave_struct.Freq(i, :), -phase_U_deg(i, :)', 'Color', Int_colors(i, :), 'LineStyle', ':', 'LineWidth', 1.5);
        plot(ISwave_struct.Freq(i, :), -phase_dQ_deg(i, :)', 'Color', Int_colors(i, :), 'Marker', '+', 'MarkerSize', 9, 'LineStyle', '-.', 'LineWidth', 1);
    end
    ax = gca;
    ax.XScale = 'log'; % for putting the scale in log
    xlim(xlim_array)
    xlabel('Frequency [Hz]');
    ylabel('Z Phase [deg]');
    legend(flipud(h), legend_flip)
    legend boxoff

%------------- END OF CODE --------------
