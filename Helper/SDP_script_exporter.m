function SDP_script_exporter(Jtr_time, prefix, varargin)
%SDP_SCRIPT_EXPORTER - Exports data of a set of Step-Dwell-Probe simulations to text files
% save the main data from a set of SDP_result structs created by
% SDP_script to text files, for easing the import with Origin (from OriginLab).
% 
% Syntax:  SDP_script_exporter(Jtr_time, prefix, SDP_results1, SDP_results2, SDP_results3)
%
% Inputs:
%   JTR_TIME - float, the time point (during the probe pulse step) that
%     should be plotted. Ideally should be enough time for the free charges
%     to reach steady state but not enough time for the ionic (which can
%     also be frozen in the SDP_script). For getting the last timepoint of
%     the applied pulse you can use Inf
%   PREFIX - char array, prefix to be used for the text files names
%   SDP_RESULTS - a struct containing the most important results of the IS simulation
%
% Example:
%   SDP_script_exporter(1e-7, 'SDP_pedot_06V', sdpsol_pedot_alt_dark, sdpsol_pedot_alt_01sun, sdpsol_pedot_alt_1sun)
%     save data from a set of simulations to text files
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also SDP_script.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

%% create header

bias_ints = zeros(length(varargin),1);

data = zeros(size(varargin{1}.Jtr, 2), 2*length(varargin));

for i = 1:length(varargin)
    assert(~isfinite(Jtr_time) || varargin{i}.t_Jtr(end) >= Jtr_time,...
        [mfilename ' - The requested time is over the pulse duration. For using always the last timepoint use Inf'])
    p1_list = find(varargin{i}.t_Jtr <= Jtr_time);
    p1(i) = p1_list(end);
    if i ~= 1
        if ~ismembertol(varargin{i-1}.t_Jtr(p1(i-1)), varargin{i}.t_Jtr(p1(i)))
            warning('Solutions are saved at different times!!!!')
            disp(varargin{i-1}.t_Jtr(p1(i-1)))
            disp(varargin{i}.t_Jtr(p1(i)))
        end
    end
    SDP_result = varargin{i};
    bias_ints(i) = SDP_result.bias_int;
    data(:,i*2-1) = SDP_result.tdwell_arr';
    data(:,i*2) = SDP_result.Jtr(p1(i),:)';
end

% round to two significant digits
legend = string(round(bias_ints, 2, 'significant'));
% add sun to numbers in legend
legend = strcat(legend, ' sun');
% replace zero in legend with dark
legend(legend=="0 sun") = "dark";

%% get header

header = strings([1,2*length(legend)]);
for i = 1:length(legend)
    header(i*2-1) = "Dwell time";
    header(i*2) = legend(i);
end

%% get measure units

units = repmat(["s", "mA/cm\+(2)"], 1,length(legend));

%% join fields

toBeSavedData = [header; units; data];

%% set NaNs to NaN

toBeSavedData = fillmissing(toBeSavedData, 'constant', "NaN");

%% save csv

fid_data = fopen([prefix '-SDP_Jtr.txt'], 'wt+');

for i = 1:size(toBeSavedData, 1)
    fprintf(fid_data, '%s\t', toBeSavedData(i, 1:end-1));
    fprintf(fid_data, '%s', toBeSavedData(i, end));
    fprintf(fid_data, '\n');
end

fclose(fid_data);

%------------- END OF CODE --------------