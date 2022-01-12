function savefig_todir(H, dir_file_name, varargin)
%SAVEFIG_TODIR - Plots a list of impedance spectroscopy (IS) results with explicit line makeup
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

valid_name = true;
try    
    java.io.File(dir_file_name).toPath;
catch
    valid_name = false;
end
if ~isempty(dir_file_name) && valid_name
    dirname = fileparts(dir_file_name);
    if ~isempty(dirname) && 7~=exist(dirname,'dir')
        mkdir(dirname)
    end
    savefig(H, [char(dir_file_name) char([varargin{:} '.fig'])], 'compact')
    saveas(H, [char(dir_file_name) char([varargin{:} '.png'])])
else
    error(['The provided path "' char(dir_file_name) '" is invalid, cannot save figure'])
end

%------------- END OF CODE --------------