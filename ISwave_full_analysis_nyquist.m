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

% create a color array with one color more than necessary
jet_matrix = jet(length(ISwave_struct.Int) + 1);
% find the yellow (which in RGB code is 1,1,0) and remove it from the
% colors list
jet_yellow_logical = ismember(jet_matrix, [1, 1, 0], 'rows');
jet_no_yellow = jet_matrix(~jet_yellow_logical, :);
Int_colors = colormap(jet_no_yellow);

% flip array and matrixes for starting from dark
legend_Int = string(flipud(ISwave_struct.Int));
re_flip = flipud(ISwave_struct.impedance_re);
im_flip = flipud(ISwave_struct.impedance_im);

% add sun to numbers in legend
legend_Int = strcat(legend_Int, ' sun');

% replace zero in legend with dark
legend_Int(legend_Int=="0 sun") = "dark";

figure('Name', 'Nyquist plot of IS at various light intensities', 'NumberTitle', 'off')
    hold off
    for i = 1:length(ISwave_struct.Int)
        plot(re_flip(i, :), im_flip(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
    end
    xlabel('Re(Z) [\Omega cm^2]');
    ylabel('-Im(Z) [\Omega cm^2]');
    legend(legend_Int)
    
%------------- END OF CODE --------------
