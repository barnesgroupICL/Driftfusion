function IS_full_analysis_impedance(IS_struct)
%IS_FULL_ANALYSIS_VSFREQUENCY - Plot Impedance Spectroscopy (IS) in a range
% of background light intensities, capacitance versus voltage oscillation
% frequency at various light bias
%
% Syntax:  IS_full_analysis_vsfrequency(IS_struct)
%
% Inputs:
%   IS_STRUCT - a struct containing the most important results of the ISstep or the ISwave simulation
%
% Example:
%   IS_full_analysis_vsfrequency(IS_struct)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also ISstep_full_exec, ISwave_full_exec.

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
if numel(unique(IS_struct.Int)) > 1
    legend_text = IS_struct.Int;
    legend_append = ' sun';
else
    legend_text = IS_struct.Vdc;
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

%% plot apparent capacitance vs frequency
figure('Name', 'IS at light intensities', 'NumberTitle', 'off');
    hold off
    for i = 1:length(legend_text)
        h(i) = plot(IS_struct.Freq(i, :), IS_struct.cap(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        if isfield(IS_struct, 'cap_idrift')
            % capacitance due to ions
            plot(IS_struct.Freq(i, :), IS_struct.cap_idrift(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--');
        end
    end
    ax = gca;
    ax.XScale = 'log'; % for putting the scale in log
    ax.YScale = 'log'; % for putting the scale in log
    xlim([min(min(IS_struct.Freq)), max(max(IS_struct.Freq))])
    xlabel('Frequency [Hz]');
    ylabel('Im(\omega^{-1} Z^{-1}) [F/cm^2]');
    legend(flipud(h), legend_flip)

%% in case of ISwave do additional graphics
if isfield(IS_struct, 'J_phase') % just ISwave have phase output

    figure('Name', 'imaginary impedance at light intensities', 'NumberTitle', 'off');
        hold off
        for i = 1:length(legend_text)
            h(i) = plot(IS_struct.Freq(i, :), IS_struct.impedance_im(i, :)',...
                'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
                'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
                'MarkerSize', 3, 'LineWidth', 1.3);
            hold on
            plot(IS_struct.Freq(i, :), IS_struct.impedance_i_im(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--');
        end
        ax = gca;
        ax.XScale = 'log'; % for putting the scale in log
        ax.YScale = 'log'; % for putting the scale in log
        xlim([min(min(IS_struct.Freq)), max(max(IS_struct.Freq))])
        xlabel('Frequency [Hz]');
        ylabel('-Im(Z) [\Omega cm^2]');
        legend(flipud(h), legend_flip)
        
    figure('Name', 'real impedance at light intensities', 'NumberTitle', 'off');
        hold off
        for i = 1:length(legend_text)
            h(i) = plot(IS_struct.Freq(i, :), IS_struct.impedance_re(i, :)',...
                'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
                'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
                'MarkerSize', 3, 'LineWidth', 1.3);
            hold on
            plot(IS_struct.Freq(i, :), IS_struct.impedance_i_re(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--');
        end
        ax = gca;
        ax.XScale = 'log'; % for putting the scale in log
        ax.YScale = 'log'; % for putting the scale in log
        xlim([min(min(IS_struct.Freq)), max(max(IS_struct.Freq))])
        xlabel('Frequency [Hz]');
        ylabel('Re(Z) [\Omega cm^2]');
        legend(flipud(h), legend_flip)
       
    figure('Name', 'abs impedance at light intensities', 'NumberTitle', 'off');
        hold off
        for i = 1:length(legend_text)
            h(i) = plot(IS_struct.Freq(i, :), IS_struct.impedance_abs(i, :)',...
                'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
                'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
                'MarkerSize', 3, 'LineWidth', 1.3);
            hold on
            plot(IS_struct.Freq(i, :), IS_struct.impedance_i_abs(i, :)', 'Color', Int_colors(i, :), 'LineStyle', '--');
        end
        ax = gca;
        ax.XScale = 'log'; % for putting the scale in log
        ax.YScale = 'log'; % for putting the scale in log
        xlim([min(min(IS_struct.Freq)), max(max(IS_struct.Freq))])
        xlabel('Frequency [Hz]');
        ylabel('Abs(Z) [\Omega cm^2]');
        legend(flipud(h), legend_flip)
end

%------------- END OF CODE --------------
