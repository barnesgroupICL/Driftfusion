function ISwave_full_analysis_nyquist(ISwave_results)
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
set(0, 'defaultAxesFontSize', 24);
% set image dimension
set(0, 'defaultfigureposition', [0, 0, 1000, 750]);
% set line thickness
set(0, 'defaultLineLineWidth', 2);

% check which was the variable being explored
if numel(unique(ISwave_results.Int)) > 1
    legend_text = ISwave_results.Int;
    legend_append = ' sun';
else
    legend_text = ISwave_results.Vdc;
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
        % find the points of the decades
        temp_logfreq = log10(ISwave_results.Freq(i,:));
        evendecade_index = mod(temp_logfreq, 2) == 0;
        
        h(i) = plot(ISwave_results.impedance_re(i, :), -ISwave_results.impedance_im(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        
        % add bigger markers for the points of the decades
        plot(ISwave_results.impedance_re(i, evendecade_index), -ISwave_results.impedance_im(i, evendecade_index)',...
            'MarkerFaceColor', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'Marker', 'x', 'LineStyle', 'none', 'MarkerSize', 12);
    end
    xlabel('Re(Z) [\Omega cm^2]');
    ylabel('-Im(Z) [\Omega cm^2]');
    legend(flipud(h(h~=0)), legend_flip)
    legend boxoff

% preallocate figures handles
h = zeros(length(legend_text), 1);
% function for normalization fitting
%semicircle = @(coeff, re) real((coeff(1).^2 - (re - coeff(2)).^2).^0.5);
figure('Name', 'Normalized Nyquist plot of IS at various light intensities', 'NumberTitle', 'off')
    hold off
    % take the last point orizontally
    %norm_array = ISwave_struct.impedance_re(:, end) / min(ISwave_struct.impedance_re(:, end));
    % take the maximum point vertically
    max_array = max(-ISwave_results.impedance_im, [], 2);
    norm_array = max_array / min(max_array);
    for i = 1:length(legend_text)
%         fit_re = ISwave_struct.impedance_re(i, end-9:end);
%         fit_im = -ISwave_struct.impedance_im(i, end-9:end);
%         try
%             fit = fitnlm(fit_re, fit_im, semicircle, [max_array(i).^2, max_array(i)]);
%             norm_array(i) = sqrt(fit.Coefficients.Estimate(1));
%         catch
%             disp([mfilename ' - Failed fit for ' legend_flip(i)])
%         end
        
        % find the points of the decades
        temp_logfreq = log10(ISwave_results.Freq(i,:));
        evendecade_index = mod(temp_logfreq, 2) == 0;
        
        h(i) = plot(ISwave_results.impedance_re(i, :) / norm_array(i), -ISwave_results.impedance_im(i, :)' / norm_array(i),...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        
        % add bigger markers for the points of the decades
        plot(ISwave_results.impedance_re(i, evendecade_index) / norm_array(i), -ISwave_results.impedance_im(i, evendecade_index)' / norm_array(i),...
            'MarkerFaceColor', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'Marker', 'x', 'LineStyle', 'none', 'MarkerSize', 12);
    end
    xlabel('Re(Z) [a.u.]');
    ylabel('-Im(Z) [a.u.]');
    % add the normalization to the legend
    legend_flip_norm = strcat(legend_flip, " / ", flipud(num2str(round(norm_array, 2, 'significant'))));
    % remove unneeded spaces from legend
    legend_flip_norm = regexprep(legend_flip_norm, ' +', ' ');
    legend(flipud(h(h~=0)), legend_flip_norm)
    legend boxoff
    
%------------- END OF CODE --------------
