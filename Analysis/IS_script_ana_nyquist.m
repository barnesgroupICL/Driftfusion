function IS_script_ana_nyquist(IS_results)
%IS_SCRIPT_ANA_NYQUIST - Plot Nyquist graph for Impedance Spectroscopy (IS)
% in a range of background light intensities or applied DC voltages, imaginary part of impedance versus real part of impedance
% Nyquist plot refers the imaginary component of the impedance to its real component.
% A normalized spectra is also plotted and the rescaling factors are
% indicated in the legend.
%
% Syntax:  IS_script_ana_nyquist(IS_results)
%
% Inputs:
%   IS_RESULTS - a struct containing the most important results of the IS simulation
%
% Example:
%   IS_script_ana_nyquist(IS_oc)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also IS_script, IS_script_ana_phase, IS_script_ana_impedance.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

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

figure('Name', 'Nyquist plot of IS at various light intensities', 'NumberTitle', 'off')
    hold off
    for i = 1:length(legend_text)
        % find the points of the even decades
        temp_logfreq = log10(IS_results.Freq(i,:));
        evendecade_index = mod(temp_logfreq, 2) == 0;
        
        h(i) = plot(IS_results.impedance_re(i, :), -IS_results.impedance_im(i, :)',...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        
        % add bigger markers for the points of the even decades
        plot(IS_results.impedance_re(i, evendecade_index), -IS_results.impedance_im(i, evendecade_index)',...
            'MarkerFaceColor', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'Marker', 'x', 'LineStyle', 'none', 'MarkerSize', 12);
    end
    xlabel('Re(Z) [\Omega cm^2]');
    ylabel('-Im(Z) [\Omega cm^2]');
    legend(flipud(h(h~=0)), legend_flip)
    legend boxoff

% preallocate figures handles
h = zeros(length(legend_text), 1);
% normalized Nyquist plot
figure('Name', 'Normalized Nyquist plot of IS at various light intensities', 'NumberTitle', 'off')
    hold off
    % take the maximum point vertically
    max_array = max(-IS_results.impedance_im, [], 2);
    norm_array = max_array / min(max_array);
    for i = 1:length(legend_text)
        
        % find the points of the even decades
        temp_logfreq = log10(IS_results.Freq(i,:));
        evendecade_index = mod(temp_logfreq, 2) == 0;
        
        h(i) = plot(IS_results.impedance_re(i, :) / norm_array(i), -IS_results.impedance_im(i, :)' / norm_array(i),...
            'Color', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'MarkerFaceColor', Int_colors(i, :), 'Marker', 's',...
            'MarkerSize', 3, 'LineWidth', 1.3);
        hold on
        
        % add bigger markers for the points of the even decades
        plot(IS_results.impedance_re(i, evendecade_index) / norm_array(i), -IS_results.impedance_im(i, evendecade_index)' / norm_array(i),...
            'MarkerFaceColor', Int_colors(i, :), 'MarkerEdgeColor', Int_colors(i, :),...
            'Marker', 'x', 'LineStyle', 'none', 'MarkerSize', 12);
    end
    xlabel('Re(Z) [\Omega cm^2]');
    ylabel('-Im(Z) [\Omega cm^2]');
    % add the normalization to the legend
    legend_flip_norm = strcat(legend_flip, " / ", flipud(num2str(round(norm_array, 2, 'significant'))));
    % remove unneeded spaces from legend
    legend_flip_norm = regexprep(legend_flip_norm, ' +', ' ');
    legend(flipud(h(h~=0)), legend_flip_norm)
    legend boxoff
    
%------------- END OF CODE --------------
