function fromISwaveResultsToTxt(ISwave_results, prefix)
% save the main data from an ISwave_struct created by ISwave_full_exec* to
% txt files, ideally easy to import with Origin (from OriginLab)

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
headerNyquist = ['Frequency', headerReIm'];

%% get measure units

unitsCap = ['Hz', repelem("F/cm\+(2)", length(legendImpedance))];
unitsNyquist = ['Hz', repelem("\g(W)·cm\+(2)", 2*length(legendImpedance))];

%% get data

% all the lines are the same for ISwave
frequencies = ISwave_results.Freq(1, :)';

cap = ISwave_results.cap';
dataCap = [frequencies, cap];

capIonic = ISwave_results.cap_idrift';
dataCapIonic = [frequencies, capIonic];

impedance_re = ISwave_results.impedance_re';
impedance_im = ISwave_results.impedance_im';
impedance = impedance_im(:, [1;1]*(1:size(impedance_im, 2)));
impedance(:, 1:2:end) = impedance_re;
dataNyquist = [frequencies, impedance];

%% join fields

toBeSavedCap = [headerFrequencyIntVdc; unitsCap; dataCap];

toBeSavedCapIonic = [headerFrequencyIntVdc; unitsCap; dataCapIonic];

toBeSavedNyquist = [headerNyquist; unitsNyquist; dataNyquist];

%% save csv

fid_cap = fopen([prefix '-cap.txt'],'wt+');
fid_capIonic = fopen([prefix '-cap_ionic.txt'],'wt+');
fid_nyquist = fopen([prefix '-nyquist.txt'],'wt+');

for i = 1:size(toBeSavedCap, 1)
    fprintf(fid_cap, '%s\t', toBeSavedCap(i, :));
    fprintf(fid_cap, '\n');
    
    fprintf(fid_capIonic, '%s\t', toBeSavedCapIonic(i, :));
    fprintf(fid_capIonic, '\n');
    
    fprintf(fid_nyquist, '%s\t', toBeSavedNyquist(i, :));
    fprintf(fid_nyquist, '\n');
end

fclose(fid_cap);
fclose(fid_capIonic);
fclose(fid_nyquist);
