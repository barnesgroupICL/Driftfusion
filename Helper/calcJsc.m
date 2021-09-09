function [EgArr, Jsc_vs_Eg] = calcJsc
% CALCJSC calcualtes the maximum theoretical short circuit current 
% as a function of band gap (Eg) using a step function absorption spectrum
% T - temperature of the device
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Start code
figson = 1;

%load AM1.5
AM15_data = xlsread('AM15.xls');
AM15_data = AM15_data';

if figson == 1
% Plot AM15 spectrum
    figure(900)
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
    figure(901)
    plot(E_ph,  AM15)
    xlabel('Energy [eV]');
    ylabel('Photon  flux density [cm-2s-1nm-1]');
    
    figure(902)
    plot(EgArr, Jsc_vs_Eg)
    xlabel('Energy [eV]');
    ylabel('J_{SC,max} [mAcm-2]')
    
    figure(903)
    plotyy(E_ph, AM15, EgArr, Jsc_vs_Eg)
    xlabel('Energy [eV]'); 
end

end