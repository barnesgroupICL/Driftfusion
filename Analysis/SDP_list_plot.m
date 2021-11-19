function SDP_list_plot(Jtr_time, dir_file_name, varargin)
%SDP_LIST_PLOT - Plots a list of Step-Dwell-Probe (SDP) results with explicit line makeup
% The plots can be saved to files if the dir_file_name option is provided.
% This is an alternative version of the anasdp.m and SDP_ana.m code.
%
% Syntax:  SDP_list_plot(Jtr_time, dir_file_name, varargin)
%
% Inputs:
%   JTR_TIME - float, the time point (during the probe pulse step) that
%     should be plotted. Ideally should be enough time for the free charges
%     to reach steady state but not enough time for the ionic (which can
%     also be frozen in the SDP_script)
%   DIR_FILE_NAME - char array, images with this prefix will be created in
%     the current directory. A path specifying the destination
%     subdirectory can also be provided. To avoid the saving, an empty char
%     array can be specified: ''.
%   VARARGIN - many arguments, the script will group the arguments in
%     groups of three:
%       the first of each group as the SDP struct with the SDP simulation
%         data, it has to contain at least the properties t_Jtr, Jtr and
%         tdwell_arr;
%       the second of each group as the legend entry;
%       the third, which has to be a {cell}, as the options to be passed to
%         the plot command (e.g. specifying the color of the line).
%
% Example:
%   SDP_list_plot(1e-7, '20210302_spiro',...
%       sdpsol_spiro_lower_1sun, 'lower 1 sun',{':r','LineWidth',3},...
%       sdpsol_spiro_1sun, '1 sun',{'--r'},...
%       sdpsol_spiro_higher_1sun, 'higher 1 sun',{'-r'})
%     plot 3 different SDP simulations, all in red but with different line
%     style and different legends. And save the graphics as files.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also SDP_script, anasdp, SDP_ana.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

for i = 1:(length(varargin)/3)
    isol = round((i*3)-2);
    p1_list = find(varargin{isol}.t_Jtr >= Jtr_time);
    p1(i) = p1_list(1);
    Jtr_time_arr(i,:) = varargin{isol}.Jtr(p1(i),:);
    if ~ismembertol(i, 1)
        if ~ismembertol(varargin{isol-3}.t_Jtr(p1(i-1)), varargin{isol}.t_Jtr(p1(i)))
            warning('Solutions are plotted at different times!!!!')
            disp(varargin{isol-3}.t_Jtr(p1(i-1)))
            disp(varargin{isol}.t_Jtr(p1(i)))
        end
    end
end

fig = figure('Name', 'SDP vs dwell time', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
hold off
ymin = Inf;
ymax = -Inf;
%legendarr = zeros(length(varargin)/3, 1);
for i = 1:(length(varargin)/3)
    isol = i*3-2;
    legendarr{i} = varargin{isol+1};
    options = varargin{isol+2};
    plot(1./varargin{isol}.tdwell_arr, Jtr_time_arr(i,:), options{:})
    ymin = min(ymin, min(Jtr_time_arr(i,ceil(end/2):end)));
    ymax = max(ymax, max(Jtr_time_arr(i,1:ceil(end/2))));
    hold on
end
xlabel('1/t_{dwell} [s-1]')
ylabel('Jtr [Acm-2]')
ax = gca;
ax.XScale = 'log'; % for putting the scale in log
%xlim([2e-8, 1e3])
range = ymax-ymin;
ylim([ymin-0.03*range, ymax+0.03*range])
legend(legendarr,'Location','southeast')
legend boxoff


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
    saveas(fig, [char(dir_file_name) char('-SDP.fig')])
    saveas(fig, [char(dir_file_name) char('-SDP.png')])
end

end

%------------- END OF CODE --------------
