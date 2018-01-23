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
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

Int_colors = colormap(jet(length(ISwave_struct.Int)));
legend_Int = string(ISwave_struct.Int);

figure('Name', 'Nyquist plot of IS at various light intensities', 'NumberTitle', 'off')
    hold off
    for i = 1:length(ISwave_struct.Int)
        plot(ISwave_struct.impedance_re(i, :), ISwave_struct.impedance_im(i, :)', 'Color', Int_colors(i, :), 'Marker', 'd');
        hold on
    end
    xlabel('Re(Z) [\Omega cm^2]');
    ylabel('Im(Z) [\Omega cm^2]');
    legend(legend_Int)
    
%------------- END OF CODE --------------
