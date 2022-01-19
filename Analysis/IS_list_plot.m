function IS_list_plot(type, dir_file_name, varargin)
%IS_LIST_PLOT - Plots a list of impedance spectroscopy (IS) results with explicit line makeup
% The plots can be saved to files if the dir_file_name option is provided
%
% Syntax:  IS_list_plot(type, dir_file_name, varargin)
%
% Inputs:
%   TYPE - char array, which quantity to plot, can be either
%     'impedance_re', 'impedance_im', 'impedance_abs', 'capacitance',
%     'phase' or 'nyquist'
%   DIR_FILE_NAME - char array, images with this prefix will be created in
%     a the current directory. A path specifying the destination
%     subdirectory can also be provided. To avoid the saving, an empty char
%     array can be specified: ''.
%   VARARGIN - many arguments, the script will group the arguments in
%     groups of three:
%       the first of each group as the IS struct with the IS simulation
%         data;
%       the second of each group as the legend entry;
%       the third, which has to be a {cell}, as the options to be passed to
%         the plot command (e.g. specifying the color of the line).
%
% Example:
%   IS_list_plot('capacitance', '20210302_spiro',...
%       IS_spiro_sc_20mV_1sun, ' 1 sun',{'--r'},...
%       IS_spiro_10x_sc_20mV_1sun, '10x 1 sun',{':r','LineWidth',3},...
%       IS_spiro_100x_sc_20mV_1sun, '100x 1 sun',{'-r'})
%     plot 3 different IS simulations, all in red but with different line
%     style and different legends. And save the graphics as files.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also IS_script, IS_script_plot_phase, IS_script_plot_impedance, IS_script_plot_nyquist.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

for i = 1:(length(varargin)/3)
    isol = i*3-2;
    legendarr{i} = varargin{isol+1};
    options = varargin{isol+2};
    switch type
        case 'impedance_re'
            fig = figure('Name', 'Real impedance Re(Z)', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
            hold off
            h(i) = plot(varargin{isol}.Freq, varargin{isol}.impedance_re, options{:});
            xlabel('Frequency [Hz]');
            ylabel('Re(Z) [\Omega cm^2]');
            ax = gca;
            ax.XScale = 'log'; % for putting the scale in log
            ax.YScale = 'log'; % for putting the scale in log
        case 'impedance_im'
            fig = figure('Name', 'Imaginary impedance -Im(Z)', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
            hold off
            h(i) = plot(varargin{isol}.Freq, -varargin{isol}.impedance_im, options{:});
            xlabel('Frequency [Hz]');
            ylabel('-Im(Z) [\Omega cm^2]');
            ax = gca;
            ax.XScale = 'log'; % for putting the scale in log
            ax.YScale = 'log'; % for putting the scale in log
        case 'impedance_abs'
            fig = figure('Name', 'Absolute impedance magnitude |Z|', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
            hold off
            h(i) = plot(varargin{isol}.Freq, varargin{isol}.impedance_abs, options{:});
            xlabel('Frequency [Hz]');
            ylabel('Abs(Z) [\Omega cm^2]');
            ax = gca;
            ax.XScale = 'log'; % for putting the scale in log
            ax.YScale = 'log'; % for putting the scale in log
        case 'capacitance'
            fig = figure('Name', 'Apparent capacitance from EIS', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
            hold off
            h(i) = plot(varargin{isol}.Freq, varargin{isol}.cap, options{:});
            hold on
            plot(varargin{isol}.Freq, -varargin{isol}.cap, options{:}, 'Marker', 'o', 'MarkerSize', 7)
            xlabel('Frequency [Hz]');
            ylabel('\omega^{-1} x Im(Z^{-1}) [F/cm^2]');
            ax = gca;
            ax.XScale = 'log'; % for putting the scale in log
            ax.YScale = 'log'; % for putting the scale in log
        case 'phase'
            fig = figure('Name', 'Phase of EIS Bode plot', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
            hold off
            phase_n_deg = rad2deg(varargin{isol}.Jtot_phase);
            h(i) = plot(varargin{isol}.Freq, -phase_n_deg, options{:});
            xlabel('Frequency [Hz]');
            ylabel('Z Phase [deg]');
            ax = gca;
            ax.XScale = 'log'; % for putting the scale in log
        case 'nyquist'
            fig = figure('Name', 'Nyquist plot of EIS', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
            hold off
            h(i) = plot(varargin{isol}.impedance_re, -varargin{isol}.impedance_im, options{:});
            xlabel('Re(Z) [\Omega cm^2]');
            ylabel('-Im(Z) [\Omega cm^2]');
        otherwise
            error([mfilename(1) ' - *type* option needed'])
    end
    %ymin = min(ymin, min());
    %ymax = max(ymax, max());
    %range = ymax-ymin;
    %ylim([ymin-0.03*range, ymax+0.03*range])
    %xlim([min(min(EA_results.Freq)), max(max(EA_results.Freq))])
    hold on
end

legend(h, legendarr, 'Location', 'northeast')
legend boxoff

savefig_todir(fig, dir_file_name, ['-IS_' type]);

end

%------------- END OF CODE --------------
