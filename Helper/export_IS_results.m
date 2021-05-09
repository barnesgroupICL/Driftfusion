function export_IS_results(IS_results, prefix)
%EXPORT_IS_RESULTS - Exports data of a set of impedance simulations to text files
% save the main data from an IS_results struct created by
% IS_script to text files, for easing the import with Origin (from OriginLab).
% 
%
% Syntax:  export_IS_results(IS_results, prefix)
%
% Inputs:
%   IS_RESULTS - a struct containing the most important results of the IS simulation
%   PREFIX - char array, prefix to be used for the text files names
%
% Example:
%   export_IS_results(IS_oc, 'IS_oc')
%     save data from a set of simulations to text files
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also IS_EA, IS_script.
% October 2017; Last revision: January 2018

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

%% create header

% check which was the variable being explored
if numel(unique(IS_results.Int)) > 1
    legend_text = IS_results.Int;
    legend_append = ' sun';
else
    legend_text = IS_results.Vdc;
    legend_append = ' Vdc';
end

% round to two significant digits
legendImpedance = round(legend_text, 2, 'significant');
% this will start from 1 sun and go to dark
legendImpedance = string(legendImpedance);
% add sun to numbers in legend
legendImpedance = strcat(legendImpedance, legend_append);
% replace zero in legend with dark
legendImpedance(legendImpedance=="0 sun") = "dark";

headerFrequencyIntVdc = ['Frequency', legendImpedance'];
headerReIm = repelem(strcat(legendImpedance, ' real'), 2);
headerReIm(2:2:end) = legendImpedance;
if iscolumn(headerReIm)
    headerNyquist = ['Frequency', headerReIm'];
else
    headerNyquist = ['Frequency', headerReIm];
end

%% get measure units

unitsCap = ['Hz', repelem("F/cm\+(2)", length(legendImpedance))];
unitsNyquist = ['Hz', repelem("\g(W)cm\+(2)", 2*length(legendImpedance))];
unitsZabs = ['Hz', repelem("\g(W)cm\+(2)", length(legendImpedance))];
unitsPhase = ['Hz', repelem("rad", length(legendImpedance))];
unitsJ = ['Hz', repelem("A/cm\+(2)", length(legendImpedance))];

%% add comments, like VOCs

commentV_DC = ["V_DC", IS_results.Vdc'];

%% get data

% all the lines are the same for IS
frequencies = IS_results.Freq(1, :)';

% capacitance
cap = IS_results.cap';
dataCap = [frequencies, cap];

% ionic capacitance
capIonic = IS_results.cap_ion_disp';
dataCapIonic = [frequencies, capIonic];

% recombination capacitance
capRec = IS_results.cap_r';
dataCapRec = [frequencies, capRec];

% accumulating current capacitance
capAcc = IS_results.cap_np_dt';
dataCapAcc = [frequencies, capAcc];

% nyquist
impedance_re = IS_results.impedance_re';
impedance_im = IS_results.impedance_im';
impedance = impedance_im(:, [1;1]*(1:size(impedance_im, 2)));
impedance(:, 1:2:end) = impedance_re;
dataNyquist = [frequencies, impedance];

% absolute impedance
impedance_abs = IS_results.impedance_abs';
dataZabs = [frequencies, impedance_abs];

% phase
Zphase = -IS_results.Jtot_phase';
dataZphase = [frequencies, Zphase];

% ionic phase
ZphaseIonic = -IS_results.ion_disp_phase';
dataZphaseIonic = [frequencies, ZphaseIonic];

% absolute ionic current amplitude
JionicAmpAbs = abs(IS_results.ion_disp_amp)';
dataJionicAmpAbs = [frequencies, JionicAmpAbs];

% absolute out of phase current amplitude
JampOutOfPhaseAbs = abs(IS_results.Jtot_amp.*sin(IS_results.Jtot_phase))';
dataJampOutOfPhaseAbs = [frequencies, JampOutOfPhaseAbs];

% absolute out of phase recombination current amplitude
JrecAmpOutOfPhaseAbs = abs(IS_results.r_amp.*sin(IS_results.r_phase))';
dataJrecAmpOutOfPhaseAbs = [frequencies, JrecAmpOutOfPhaseAbs];

%% join fields

toBeSavedCap = [headerFrequencyIntVdc; unitsCap; dataCap];

toBeSavedCapIonic = [headerFrequencyIntVdc; unitsCap; dataCapIonic];

toBeSavedCapRec = [headerFrequencyIntVdc; unitsCap; dataCapRec];

toBeSavedCapAcc = [headerFrequencyIntVdc; unitsCap; dataCapAcc];

toBeSavedNyquist = [headerNyquist; unitsNyquist; dataNyquist];

toBeSavedZabs = [headerFrequencyIntVdc; unitsZabs; dataZabs];

toBeSavedZphase = [headerFrequencyIntVdc; unitsPhase; dataZphase];

toBeSavedZphaseIonic = [headerFrequencyIntVdc; unitsPhase; dataZphaseIonic];

toBeSavedJionicAmpAbs = [headerFrequencyIntVdc; unitsJ; commentV_DC; dataJionicAmpAbs];

toBeSavedJampOutOfPhaseAbs = [headerFrequencyIntVdc; unitsJ; commentV_DC; dataJampOutOfPhaseAbs];

toBeSavedJrecAmpOutOfPhaseAbs = [headerFrequencyIntVdc; unitsJ; commentV_DC; dataJrecAmpOutOfPhaseAbs];

%% set NaNs to NaN

toBeSavedCap = fillmissing(toBeSavedCap, 'constant', "NaN");
toBeSavedCapIonic = fillmissing(toBeSavedCapIonic, 'constant', "NaN");
toBeSavedCapRec = fillmissing(toBeSavedCapRec, 'constant', "NaN");
toBeSavedCapAcc = fillmissing(toBeSavedCapAcc, 'constant', "NaN");
toBeSavedNyquist = fillmissing(toBeSavedNyquist, 'constant', "NaN");
toBeSavedZabs = fillmissing(toBeSavedZabs, 'constant', "NaN");
toBeSavedZphase = fillmissing(toBeSavedZphase, 'constant', "NaN");
toBeSavedZphaseIonic = fillmissing(toBeSavedZphaseIonic, 'constant', "NaN");
toBeSavedJionicAmpAbs = fillmissing(toBeSavedJionicAmpAbs, 'constant', "NaN");
toBeSavedJampOutOfPhaseAbs = fillmissing(toBeSavedJampOutOfPhaseAbs, 'constant', "NaN");
toBeSavedJrecAmpOutOfPhaseAbs = fillmissing(toBeSavedJrecAmpOutOfPhaseAbs, 'constant', "NaN");

%% save csv

fid_cap = fopen([prefix '-cap.txt'], 'wt+');
fid_capIonic = fopen([prefix '-capIonic.txt'], 'wt+');
fid_capRec = fopen([prefix '-capRecombination.txt'], 'wt+');
fid_capAcc = fopen([prefix '-capAccumulating.txt'], 'wt+');
fid_nyquist = fopen([prefix '-nyquist.txt'], 'wt+');
fid_Zabs = fopen([prefix '-Zabs.txt'], 'wt+');
fid_phase = fopen([prefix '-Zphase.txt'], 'wt+');
fid_phaseIonic = fopen([prefix '-ZphaseIonic.txt'], 'wt+');
fid_JionicAmpAbs = fopen([prefix '-JionicAmpAbs.txt'], 'wt+');
fid_JampOutOfPhaseAbs = fopen([prefix '-JampOutOfPhaseAbs.txt'], 'wt+');
fid_JrecAmpOutOfPhaseAbs = fopen([prefix '-JrecAmpOutOfPhaseAbs.txt'], 'wt+');

for i = 1:size(toBeSavedCap, 1)
    fprintf(fid_cap, '%s\t', toBeSavedCap(i, 1:end-1));
    fprintf(fid_cap, '%s', toBeSavedCap(i, end));
    fprintf(fid_cap, '\n');
    
    fprintf(fid_capIonic, '%s\t', toBeSavedCapIonic(i, 1:end-1));
    fprintf(fid_capIonic, '%s', toBeSavedCapIonic(i, end));
    fprintf(fid_capIonic, '\n');
    
    fprintf(fid_capRec, '%s\t', toBeSavedCapRec(i, 1:end-1));
    fprintf(fid_capRec, '%s', toBeSavedCapRec(i, end));
    fprintf(fid_capRec, '\n');
    
    fprintf(fid_capAcc, '%s\t', toBeSavedCapAcc(i, 1:end-1));
    fprintf(fid_capAcc, '%s', toBeSavedCapAcc(i, end));
    fprintf(fid_capAcc, '\n');
    
    fprintf(fid_nyquist, '%s\t', toBeSavedNyquist(i, 1:end-1));
    fprintf(fid_nyquist, '%s', toBeSavedNyquist(i, end));
    fprintf(fid_nyquist, '\n');
    
    fprintf(fid_Zabs, '%s\t', toBeSavedZabs(i, 1:end-1));
    fprintf(fid_Zabs, '%s', toBeSavedZabs(i, end));
    fprintf(fid_Zabs, '\n');
    
    fprintf(fid_phase, '%s\t', toBeSavedZphase(i, 1:end-1));
    fprintf(fid_phase, '%s', toBeSavedZphase(i, end));
    fprintf(fid_phase, '\n');
    
    fprintf(fid_phaseIonic, '%s\t', toBeSavedZphaseIonic(i, 1:end-1));
    fprintf(fid_phaseIonic, '%s', toBeSavedZphaseIonic(i, end));
    fprintf(fid_phaseIonic, '\n');
    
    fprintf(fid_JionicAmpAbs, '%s\t', toBeSavedJionicAmpAbs(i, 1:end-1));
    fprintf(fid_JionicAmpAbs, '%s', toBeSavedJionicAmpAbs(i, end));
    fprintf(fid_JionicAmpAbs, '\n');
    
    fprintf(fid_JampOutOfPhaseAbs, '%s\t', toBeSavedJampOutOfPhaseAbs(i, 1:end-1));
    fprintf(fid_JampOutOfPhaseAbs, '%s', toBeSavedJampOutOfPhaseAbs(i, end));
    fprintf(fid_JampOutOfPhaseAbs, '\n');
    
    fprintf(fid_JrecAmpOutOfPhaseAbs, '%s\t', toBeSavedJrecAmpOutOfPhaseAbs(i, 1:end-1));
    fprintf(fid_JrecAmpOutOfPhaseAbs, '%s', toBeSavedJrecAmpOutOfPhaseAbs(i, end));
    fprintf(fid_JrecAmpOutOfPhaseAbs, '\n');
end

fclose(fid_cap);
fclose(fid_capIonic);
fclose(fid_capRec);
fclose(fid_capAcc);
fclose(fid_nyquist);
fclose(fid_Zabs);
fclose(fid_phase);
fclose(fid_phaseIonic);
fclose(fid_JionicAmpAbs);
fclose(fid_JampOutOfPhaseAbs);
fclose(fid_JrecAmpOutOfPhaseAbs);

%------------- END OF CODE --------------

