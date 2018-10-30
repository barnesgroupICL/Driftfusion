function fromISwaveResultsToTxt(ISwave_results, prefix)
% save the main data from an ISwave_results struct created by ISwave_full_exec* to txt files, ideally easy to import with Origin (from OriginLab)
% example:
%   fromISwaveResultsToTxt(ISwave_oc, 'pinParams_ISwave_oc')

%% create header

% check which was the variable being explored
if numel(unique(ISwave_results.Int)) > 1
    legend_text = ISwave_results.Int;
    legend_append = ' sun';
else
    legend_text = ISwave_results.Vdc;
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

commentV_DC = ["V_DC", ISwave_results.Vdc'];

%% get data

% all the lines are the same for ISwave
frequencies = ISwave_results.Freq(1, :)';

% capacitance
cap = ISwave_results.cap';
dataCap = [frequencies, cap];

% ionic capacitance
capIonic = ISwave_results.cap_idrift';
dataCapIonic = [frequencies, capIonic];

% recombination capacitance
capRec = ISwave_results.cap_U';
dataCapRec = [frequencies, capRec];

% accumulating current capacitance
capAcc = ISwave_results.cap_dQ';
dataCapAcc = [frequencies, capAcc];

% nyquist
impedance_re = ISwave_results.impedance_re';
impedance_im = ISwave_results.impedance_im';
impedance = impedance_im(:, [1;1]*(1:size(impedance_im, 2)));
impedance(:, 1:2:end) = impedance_re;
dataNyquist = [frequencies, impedance];

% absolute impedance
impedance_abs = ISwave_results.impedance_abs';
dataZabs = [frequencies, impedance_abs];

% phase
Zphase = -ISwave_results.J_phase';
dataZphase = [frequencies, Zphase];

% ionic phase
ZphaseIonic = -ISwave_results.J_i_phase';
dataZphaseIonic = [frequencies, ZphaseIonic];

% absolute ionic current amplitude
JionicAmpAbs = abs(ISwave_results.J_i_amp)';
dataJionicAmpAbs = [frequencies, JionicAmpAbs];

% absolute out of phase current amplitude
JampOutOfPhaseAbs = abs(ISwave_results.J_amp.*sin(ISwave_results.J_phase))';
dataJampOutOfPhaseAbs = [frequencies, JampOutOfPhaseAbs];

% absolute out of phase recombination current amplitude
JrecAmpOutOfPhaseAbs = abs(ISwave_results.J_U_amp.*sin(ISwave_results.J_U_phase))';
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
