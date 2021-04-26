function export_SDP_results(prefix, varargin)
%EXPORT_IS_RESULTS - Exports data of a set of Step-Dwell-Probe simulations to text files
% save the main data from a set of SDP_result structs created by
% doSDP_alt to text files, for easing the import with Origin (from OriginLab).
% 
%
% Syntax:  export_SDP_results(prefix, SDP_results1, SDP_results2, SDP_results3)
%
% Inputs:
%   PREFIX - char array, prefix to be used for the text files names
%   SDP_RESULTS - a struct containing the most important results of the IS simulation
%
% Example:
%   export_SDP_results('SDP_pedot_06V', sdpsol_pedot_alt_dark, sdpsol_pedot_alt_01sun, sdpsol_pedot_alt_1sun)
%     save data from a set of simulations to text files
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also doSDP_alt.

% Author: Dr. Ilario Gelmetti
% Institute of Chemical Research of Catalonia (ICIQ)
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% April 2021

%------------- BEGIN CODE --------------

%% create header

bias_ints = zeros(length(varargin),1);

data = zeros(size(varargin{1}.Jtr, 2), 2*length(varargin));

for i = 1:length(varargin)
    SDP_result = varargin{i};
    bias_ints(i) = SDP_result.bias_int;
    data(:,i*2-1) = SDP_result.tdwell_arr';
    data(:,i*2) = SDP_result.Jtr(end,:)';
end

% round to two significant digits
legendSDP = string(round(bias_ints, 2, 'significant'));
% add sun to numbers in legend
legendSDP = strcat(legendSDP, ' sun');
% replace zero in legend with dark
legendSDP(legendSDP=="0 sun") = "dark";

%% get header

header = strings([1,2*length(legendSDP)]);
for i = 1:length(legendSDP)
    header(i*2-1) = "Dwell time";
    header(i*2) = legendSDP(i);
end

%% get measure units

units = repmat(["s", "mA/cm\+(2)"], 1,length(legendSDP));

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