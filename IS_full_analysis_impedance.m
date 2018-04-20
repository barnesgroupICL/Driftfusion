function IS_full_analysis_impedance(IS_results)
%IS_FULL_ANALYSIS_VSFREQUENCY - Plot Impedance Spectroscopy (IS) in a range
% of background light intensities, capacitance versus voltage oscillation
% frequency at various light bias
%
% Syntax:  IS_full_analysis_impedance(IS_struct)
%
% Inputs:
%   IS_STRUCT - a struct containing the most important results of the ISstep or the ISwave simulation
%
% Example:
%   IS_full_analysis_impedance(IS_struct)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also ISstep_full_exec, ISwave_full_exec, ISwave_full_analysis_phase, ISwave_full_analysis_nyquist.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

% increase graphics font size
set(0, 'defaultAxesFontSize', 24);
% set image dimension
set(0, 'defaultfigureposition', [0, 0, 1000, 750]);
% set line thickness
set(0, 'defaultLineLineWidth', 2);

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

%% plot apparent capacitance vs frequency
figure('Name', 'IS at light intensities. Dashed: ionic; dotted: recombination; dashdotted: stored charge variation', 'NumberTitle', 'off');
    hold off
    for i = 1:length(legend_text)
        h(i) = plot(IS_results.Freq(i, :), IS_results.cap(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        if isfield(IS_results, 'cap_idrift')
            % capacitance due to ions
            plot(IS_results.Freq(i, :), IS_results.cap_idrift(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--', 'LineWidth', 1);
            % capacitance calculated from the recombination amount
            plot(IS_results.Freq(i, :), IS_results.cap_U(i, :)', 'Color', Int_colors(i, :), 'LineStyle', ':', 'LineWidth', 1.5);
        end
        plot(IS_results.Freq(i, :), IS_results.cap_dQ(i, :)', 'Color', Int_colors(i, :), 'Marker', '+', 'MarkerSize', 9, 'LineStyle', '-.', 'LineWidth', 1);
    end
    ax = gca;
    ax.XScale = 'log'; % for putting the scale in log
    ax.YScale = 'log'; % for putting the scale in log
    xlim(xlim_array)
    xlabel('Frequency [Hz]');
    ylabel('\omega^{-1}·Im(Z^{-1}) [F/cm^2]');
    legend(flipud(h), legend_flip)
    legend boxoff

%% in case of ISwave do additional graphics
if isfield(IS_results, 'J_phase') % just ISwave have phase output

    figure('Name', 'imaginary impedance at light intensities. Dashed: ionic; dotted: recombination; dashdotted: stored charge variation', 'NumberTitle', 'off');
        hold off
        for i = 1:length(legend_text)
            h(i) = plot(IS_results.Freq(i, :), -IS_results.impedance_im(i, :)',...
                'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
                'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
                'MarkerSize', 3, 'LineWidth', 1.3);
            hold on
            plot(IS_results.Freq(i, :), -IS_results.impedance_i_im(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--', 'LineWidth', 1);
            plot(IS_results.Freq(i, :), -IS_results.impedance_U_im(i, :)', 'Color', Int_colors(i, :), 'LineStyle', ':', 'LineWidth', 1.5);
            plot(IS_results.Freq(i, :), -IS_results.impedance_dQ_im(i, :)', 'Color', Int_colors(i, :), 'Marker', '+', 'MarkerSize', 9, 'LineStyle', '-.', 'LineWidth', 1);
        end
        ax = gca;
        ax.XScale = 'log'; % for putting the scale in log
        ax.YScale = 'log'; % for putting the scale in log
        xlim(xlim_array)
        xlabel('Frequency [Hz]');
        ylabel('-Im(Z) [\Omega cm^2]');
        legend(flipud(h), legend_flip)
        legend boxoff
        
    figure('Name', 'real impedance at light intensities. Dashed: ionic; dotted: recombination; dashdotted: stored charge variation', 'NumberTitle', 'off');
        hold off
        for i = 1:length(legend_text)
            h(i) = plot(IS_results.Freq(i, :), IS_results.impedance_re(i, :)',...
                'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
                'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
                'MarkerSize', 3, 'LineWidth', 1.3);
            hold on
            plot(IS_results.Freq(i, :), IS_results.impedance_i_re(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--', 'LineWidth', 1);
            plot(IS_results.Freq(i, :), IS_results.impedance_U_re(i, :)', 'Color', Int_colors(i, :), 'LineStyle', ':', 'LineWidth', 1.5);
            plot(IS_results.Freq(i, :), IS_results.impedance_dQ_re(i, :)', 'Color', Int_colors(i, :), 'Marker', '+', 'MarkerSize', 9, 'LineStyle', '-.', 'LineWidth', 1);
        end
        ax = gca;
        ax.XScale = 'log'; % for putting the scale in log
        ax.YScale = 'log'; % for putting the scale in log
        xlim(xlim_array)
        xlabel('Frequency [Hz]');
        ylabel('Re(Z) [\Omega cm^2]');
        legend(flipud(h), legend_flip)
        legend boxoff
       
    figure('Name', 'abs impedance at light intensities. Dashed: ionic; dotted: recombination; dashdotted: stored charge variation', 'NumberTitle', 'off');
        hold off
        for i = 1:length(legend_text)
            h(i) = plot(IS_results.Freq(i, :), IS_results.impedance_abs(i, :)',...
                'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
                'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
                'MarkerSize', 3, 'LineWidth', 1.3);
            hold on
            plot(IS_results.Freq(i, :), IS_results.impedance_i_abs(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--', 'LineWidth', 1);
            plot(IS_results.Freq(i, :), IS_results.impedance_U_abs(i, :)', 'Color', Int_colors(i, :), 'LineStyle', ':', 'LineWidth', 1.5);
            plot(IS_results.Freq(i, :), IS_results.impedance_dQ_abs(i, :)', 'Color', Int_colors(i, :), 'Marker', '+', 'MarkerSize', 9, 'LineStyle', '-.', 'LineWidth', 1);
        end
        ax = gca;
        ax.XScale = 'log'; % for putting the scale in log
        ax.YScale = 'log'; % for putting the scale in log
        xlim(xlim_array)
        xlabel('Frequency [Hz]');
        ylabel('Abs(Z) [\Omega cm^2]');
        legend(flipud(h), legend_flip)
        legend boxoff
end

%------------- END OF CODE --------------
