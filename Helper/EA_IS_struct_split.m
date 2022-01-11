function EA_IS_struct_split(in)
%EA_IS_STRUCT_SPLIT - Splits EA_script and IS_script solutions into one solution for each illumination/voltage bias
% This can be used for then plotting the individual conditions using
% EA_list_plot or IS_list_plot which do not accept complex structures with
% more than one illumination/voltage bias condition
% 
% Syntax:  EA_IS_struct_split(in)
%
% Inputs:
%   IN - a structure containing the results of a simulation performed with
%   EA_script or IS_script
%
% Outputs:
%
% Example:
%   EA_IS_struct_split(IS_sc_20mV)
%     saves the content of IS_sc_20mV in the base Matlab workspace but
%     divided in smaller structures like IS_sc_20mV_dark for the solution
%     in dark and IS_sc_20mV_sun0_1 for the solution at 0.1 sun
%     illumination
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also df.
%
%% LICENSE
% Copyright (C) 2022  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%------------- BEGIN CODE --------------

if numel(unique(in.Int)) > 1
    params = in.Int;
    scanning_illumination = true;
else
    params = in.Vdc;
    scanning_illumination = false;
end

common_content.periods = in.periods;
common_content.tpoints = in.tpoints;
common_content.deltaV = in.deltaV;
common_content.sol_name = in.sol_name;
suffixes = "";

for i = 1:length(params)
    if scanning_illumination
        if params(i)
            suffixes(i) = matlab.lang.makeValidName(['sun' num2str(params(i))]);
        else
            suffixes(i) = 'dark';
        end
    else
        suffixes(i) = matlab.lang.makeValidName(['Vapp' num2str(params(i))]);
    end
    command_generate_structure = strcat(suffixes(i), ' = common_content;');
    eval(command_generate_structure);
end

names = fieldnames(in);

for i = 1:length(names)
    command_populate_content = strcat('content = in.', names(i), ';');
    eval(command_populate_content{1});
    if size(content,1) >= 2
        disp(names(i))
        for j = 1:length(suffixes)
            command_fill_structure = strcat(suffixes(j), '.', names(i), ' = content(', num2str(j), ',:);');
            eval(command_fill_structure{1});
        end
    end
end


for i = 1:length(suffixes)
    command_save_structure = strcat("assignin('base', [inputname(1) '_", suffixes(i), "'], ", suffixes(i), ");");
    eval(command_save_structure);
end

%------------- END OF CODE --------------