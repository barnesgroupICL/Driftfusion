function EA_script_exporter(prefix, EA_results)
%EXPORT_EA_RESULTS - Exports data of a set of Electro Absorbance simulations to text files
% save the main data from an EA_results struct created by
% EA_script to text files, for easing the import with Origin (from OriginLab).
% 
%
% Syntax:  export_EA_results(prefix, EA_results)
%
% Inputs:
%   PREFIX - char array, prefix to be used for the text files names
%   EA_RESULTS - a struct containing the most important results of the EA simulation
%
% Example:
%   export_EA_results('EA_oc', EA_oc)
%     save data from a set of simulations to text files
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also doIS_EA, EA_script.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------


%                       sol_name: 'soleq_sp_p200_spD70mV_tD70mV_pepp32_spepp12_cmu01_st_Int_0'
%                            Vdc: [3×1 double]
%                        periods: 20
%                           Freq: [3×17 double]
%                        tpoints: 1601
%                           tmax: [3×17 double]
%                            Int: [3×1 double]
%                         deltaV: 0.8000
%                      sun_index: 3
%                 AC_ExDC_E_bias: [3×17 double]
%                  AC_ExDC_E_amp: [3×17 double]
%                AC_ExDC_E_phase: [3×17 double]
%               AC_ExDC_E_i_bias: [3×17 double]
%                AC_ExDC_E_i_amp: [3×17 double]
%              AC_ExDC_E_i_phase: [3×17 double]
%                AC_Efield2_bias: [3×17 double]
%                 AC_Efield2_amp: [3×17 double]
%               AC_Efield2_phase: [3×17 double]
%              AC_Efield2_i_bias: [3×17 double]
%               AC_Efield2_i_amp: [3×17 double]
%             AC_Efield2_i_phase: [3×17 double]
%     AC_Efield_amp_squared_mean: [3×17 double]


%% create header

% check which was the variable being explored
if numel(unique(EA_results.Int)) > 1
    legend_text = EA_results.Int;
    legend_append = ' sun';
else
    legend_text = EA_results.Vdc;
    legend_append = ' Vdc';
end

% round to two significant digits
legend = round(legend_text, 2, 'significant');
% this will start from 1 sun and go to dark
legend = string(legend);
% add sun to numbers in legend
legend = strcat(legend, legend_append);
% replace zero in legend with dark
legend(legend=="0 sun") = "dark";
% repeat three times each element (e.g. 1 sun dark 1 sun dark 1 sun dark)
legend_expanded = repmat(legend', 1, 3);
type = [" bias", " amp", " phase"];
type = repelem(type, length(legend));
% add type to each entry of the legend
legend = strcat(legend_expanded, type);
header = ['Frequency', legend];

%% get measure units

units = ['Hz', repelem(["V\+(2)/cm\+(2)", "V\+(2)/cm\+(2)", "rad"], round(length(legend)/3))];

%% get data

% all the lines are the same for EA
frequencies = EA_results.Freq(1, :)';

% first harmonic EA
H1 = [EA_results.AC_ExDC_E_bias',EA_results.AC_ExDC_E_amp',EA_results.AC_ExDC_E_phase'];
dataH1 = [frequencies, H1];

% second harmonic EA
H2 = [EA_results.AC_Efield2_bias',EA_results.AC_Efield2_amp',EA_results.AC_Efield2_phase'];
dataH2 = [frequencies, H2];

%% join fields
toBeSavedH1 = [header; units; dataH1];
toBeSavedH2 = [header; units; dataH2];

%% set NaNs to NaN

toBeSavedH1 = fillmissing(toBeSavedH1, 'constant', "NaN");
toBeSavedH2 = fillmissing(toBeSavedH2, 'constant', "NaN");

%% save csv

fid_H1 = fopen([prefix '-EA_H1.txt'], 'wt+');
fid_H2 = fopen([prefix '-EA_H2.txt'], 'wt+');

for i = 1:size(toBeSavedH1, 1)
    fprintf(fid_H1, '%s\t', toBeSavedH1(i, 1:end-1));
    fprintf(fid_H1, '%s', toBeSavedH1(i, end));
    fprintf(fid_H1, '\n');
end

for i = 1:size(toBeSavedH2, 1)
    fprintf(fid_H2, '%s\t', toBeSavedH2(i, 1:end-1));
    fprintf(fid_H2, '%s', toBeSavedH2(i, end));
    fprintf(fid_H2, '\n');
end

fclose(fid_H1);
fclose(fid_H2);

%------------- END OF CODE --------------

