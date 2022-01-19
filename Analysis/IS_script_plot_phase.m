function IS_script_plot_phase(IS_results)
%IS_SCRIPT_PLOT_PHASE - Represents Bode plots of phase from impedance spectroscopy
% The phase of the current oscillation with regards to the applied voltage
% is plotted as obtained from simulations made with
% IS_script or IS_script_nonparallel
%
% Syntax:  IS_script_plot_phase(IS_results)
%
% Inputs:
%   IS_RESULTS - a struct containing the most important results of the IS simulation
%
% Example:
%   IS_script_plot_phase(IS_oc)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also IS_script, IS_script_ana_nyquist, IS_script_ana_impedance, IS_list_plot.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

% check which was the variable being explored
if numel(unique(IS_results.Int)) > 1
    legend_text = IS_results.Int;
    legend_append = ' sun';
else
    legend_text = IS_results.Vdc;
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
xlim_array = [min(min(IS_results.Freq))*0.99, max(max(IS_results.Freq)*1.01)];

%% do IS phase plot

% we want to plot phase of impedance, which is -phase(I vs V)
phase_n_deg = rad2deg(IS_results.Jtot_phase);
phase_i_deg = rad2deg(IS_results.ion_disp_phase);
phase_U_deg = rad2deg(IS_results.r_phase);
phase_dQ_deg = rad2deg(IS_results.np_dt_phase);

figure('Name', 'Phase of EIS Bode plot. Dashed: ionic; dotted: recombination; dashdotted: stored charge variation', 'NumberTitle', 'off')
    hold off
    for i = 1:length(legend_text)
        h(i) = plot(IS_results.Freq(i, :), -phase_n_deg(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        % plot ionic displacement current phase
        plot(IS_results.Freq(i, :), -phase_i_deg(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--', 'LineWidth', 1);
        % plot recombination current phase
        plot(IS_results.Freq(i, :), -phase_U_deg(i, :)', 'Color', Int_colors(i, :), 'LineStyle', ':', 'LineWidth', 1.5);
        % plot accumulating current phase
        plot(IS_results.Freq(i, :), -phase_dQ_deg(i, :)', 'Color', Int_colors(i, :), 'Marker', '+', 'MarkerSize', 9, 'LineStyle', '-.', 'LineWidth', 1);
    end
    ax = gca;
    ax.XScale = 'log'; % for putting the scale in log
    xlim(xlim_array)
    xlabel('Frequency [Hz]');
    ylabel('Z Phase [deg]');
    legend(flipud(h), legend_flip)
    legend boxoff

%------------- END OF CODE --------------
