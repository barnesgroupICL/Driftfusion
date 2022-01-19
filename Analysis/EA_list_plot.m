function EA_list_plot(type, dir_file_name, varargin)
%EA_LIST_PLOT - Plots a list of electroabsorbance (EA, Stark spectroscopy) results with explicit line makeup
% The plots can be saved to files if the dir_file_name option is provided
%
% Syntax:  EA_list_plot(type, dir_file_name, varargin)
%
% Inputs:
%   TYPE - char array, which quantity to plot, can be either
%     '1h', '2h' or 'phase'
%   DIR_FILE_NAME - char array, images with this prefix will be created in
%     the current directory. A path specifying the destination
%     subdirectory can also be provided. To avoid the saving, an empty char
%     array can be specified: ''.
%   VARARGIN - many arguments, the script will group the arguments in
%     groups of three:
%       the first of each group as the EA struct with the IS simulation
%         data;
%       the second of each group as the legend entry;
%       the third, which has to be a {cell}, as the options to be passed to
%         the plot command (e.g. specifying the color of the line).
%
% Example:
%   EA_list_plot('2h', '20210409_spiro',...
%       EA_spiro_10000x_sc_800mV_1sun, '10000x 1 sun',{':r','LineWidth',3},...
%       EA_spiro_1000x_sc_800mV_1sun, '1000x 1 sun',{'-r'},...
%       EA_spiro_sc_800mV_1sun, ' 1 sun',{'--r'})
%     plot 3 different EA simulations, all in red but with different line
%     style and different legends. And save the graphics as files.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also EA_script, EA_script_plot_Efield, EA_script_plot_phase.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

switch type
    case '1h'
        fig = figure('Name', 'Amplitude of EA first harmonic E_{AC}*E_{DC}', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
    case '2h'
        fig = figure('Name', 'Amplitude of EA second harmonic E_{AC}^2', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
    case 'phase'
        fig = figure('Name', 'Phase of EA first harmonic', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
    otherwise
        error([mfilename(1) ' - *type* option needed'])
end

hold off

for i = 1:(length(varargin)/3)
    isol = i*3-2;
    legendarr{i} = varargin{isol+1};
    options = varargin{isol+2};
    switch type
        case '1h'
            plot(varargin{isol}.Freq, varargin{isol}.AC_ExDC_E_amp, options{:})
        case '2h'
            plot(varargin{isol}.Freq, varargin{isol}.AC_Efield2_amp, options{:})
        case 'phase'
            phase_n_deg = rad2deg(wrapTo2Pi(varargin{isol}.AC_ExDC_E_phase));
            plot(varargin{isol}.Freq, phase_n_deg, options{:})
    end
    hold on
end

ax = gca;
ax.XScale = 'log'; % for putting the scale in log
ax.YScale = 'log'; % for putting the scale in log
xlabel('Frequency [Hz]');
ylabel('Abs(E_{AC}^2) [V^2/cm^2]');

legend(legendarr, 'Location', 'southeast')
legend boxoff

savefig_todir(fig, dir_file_name, ['-EA_' type]);

end

%------------- END OF CODE --------------
