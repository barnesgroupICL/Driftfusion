function EA_script_plot_phase(EA_results, savefig_dir)
%EA_SCRIPT_PLOT_PHASE - Plot phase Bode plots from ElectroAbsorbance (EA)
% in a range of background light intensities or applied DC voltages
% This script is very similar to IS_script_ana_phase
% Plots the first and second harmonic phase shift with regards to the
% applied voltage.
%
% Syntax:  EA_script_plot_phase(EA_results)
%
% Inputs:
%   EA_RESULTS - a struct containing the most important results of the EA simulation
%   SAVEFIG_DIR - optional string, if provided the figures will be saved in
%     at that path
%
% Example:
%   EA_script_ana_phase(EA_results)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also EA_script, IS_script_ana_phase, EA_script_ana_Efield.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

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

if nargin > 1
    [~, legend_prepend, ~] = fileparts(savefig_dir);
    legend_prepend = regexprep(legend_prepend,'_',' ');
    % ridiculously enough I didn't find a way to add a space separator, so
    % I use a colon
    legend_flip = strcat(legend_prepend, ':', legend_flip);
end

% in case just a single frequency was simulated, min and max are
% coincident, so they are modified by a small percentage
xlim_array = [min(min(EA_results.Freq))*0.99, max(max(EA_results.Freq)*1.01)];

%% do EA first and second harmonic phase plots

phase_n_deg = rad2deg(wrapTo2Pi(EA_results.AC_ExDC_E_phase));
phase_i_deg = rad2deg(wrapTo2Pi(EA_results.AC_ExDC_E_i_phase));
fig_1h = figure('Name', 'Phase plot of EA first harmonic', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
    hold off
    for i = 1:length(legend_text)
        h(i) = plot(EA_results.Freq(i, :), phase_n_deg(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        plot(EA_results.Freq(i, :), phase_i_deg(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--');
    end
    ax = gca;
    ax.XScale = 'log'; % for putting the scale in log
    xlim(xlim_array)
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    legend(flipud(h), legend_flip)
    legend boxoff

phase_2h_deg = rad2deg(wrapToPi(EA_results.AC_Efield2_phase));
phase_2h_i_deg = rad2deg(wrapTo2Pi(EA_results.AC_Efield2_i_phase));
fig_2h = figure('Name', 'Phase plot of EA second harmonic', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
    hold off
    for i = 1:length(legend_text)
        h(i) = plot(EA_results.Freq(i, :), phase_2h_deg(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        plot(EA_results.Freq(i, :), phase_2h_i_deg(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--');
    end
    ax = gca;
    ax.XScale = 'log'; % for putting the scale in log
    xlim(xlim_array)
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    legend(flipud(h), legend_flip)
    legend boxoff

if nargin > 1
    if 7~=exist(savefig_dir,'dir')
        mkdir(savefig_dir)
    end
    pathname = [char(savefig_dir) filesep char(inputname(1))];
    saveas(fig_1h, [pathname '-phase-1h.fig'])
    saveas(fig_2h, [pathname '-phase-2h.fig'])
    saveas(fig_1h, [pathname '-phase-1h.png'])
    saveas(fig_2h, [pathname '-phase-2h.png'])
end

%------------- END OF CODE --------------
