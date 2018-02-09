function ISwave_full_analysis_nyquist(ISwave_struct)
%ISSTEP_FULL_ANALYSIS_NYQUIST - Plot Nyquist graph for Impedance Spectroscopy (ISwave) in a range
% of background light intensities, imaginary part of capacitance versus
% real part of capacitance
%
% Syntax:  ISwave_full_analysis_nyquist(ISwave_struct)
%
% Inputs:
%   ISWAVE_STRUCT - a struct containing the most important results of the ISwave simulation
%
% Example:
%   ISwave_full_analysis_nyquist(ISwave_struct)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also ISwave_full_exec.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

% increase graphics font size
set(0, 'defaultAxesFontSize', 30)

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

figure('Name', 'Nyquist plot of IS at various light intensities', 'NumberTitle', 'off')
    hold off
    for i = 1:length(legend_text)
        h(i) = plot(ISwave_struct.impedance_re(i, :), ISwave_struct.impedance_im(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
    end
    xlabel('Re(Z) [\Omega cm^2]');
    ylabel('-Im(Z) [\Omega cm^2]');
    legend(flipud(h), legend_flip)
    
%------------- END OF CODE --------------
