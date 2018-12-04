function [EgArr, Jsc_vs_Eg] = calcJsc
% CALCJSC calcualtes the maximum theoretical short circuit current 
% as a function of band gap (Eg) using a step function absorption spectrum
% T - temperature of the device

figson = 1;

%load AM1.5
AM15_data = xlsread('AM15.xls');
AM15_data = AM15_data';

if figson == 1
% Plot AM15 spectrum
    figure(600)
    plot(AM15_data(1, :), AM15_data(2, :));
    xlabel('Wavelength [nm]');
    ylabel('Power density [mWcm-2nm-1]');
    
end

% Convert AM1.5 to correct units
E_ph = 1239.8./AM15_data(1, :);            % Photon energy
AM15 = AM15_data(2, :);

E_ph = fliplr(E_ph);
AM15 = fliplr(AM15);

EgArr = 0.5:0.01:4;

for i=1:length(EgArr)
    Eg_temp = EgArr(i);
    p = find(E_ph < Eg_temp);
    p = p(end);
    
    AM15_temp = AM15(p:end);
    nm = 0.5:0.5:length(AM15_temp)*0.5;
    Jsc_vs_Eg(i) = trapz(nm, AM15_temp);
end

if figson ==1
    
    figure(601)
    plot(E_ph,  AM15)
    xlabel('Energy [eV]');
    ylabel('Photon  flux density [cm-2s-1nm-1]');
    
    figure(602)
    plot(EgArr, Jsc_vs_Eg)
    xlabel('Energy [eV]');
    ylabel('J_{SC,max} [mAcm-2]')
    
    figure(603)
    plotyy(E_ph, AM15, EgArr, Jsc_vs_Eg)
    xlabel('Energy [eV]');
    
end

end