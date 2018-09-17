function [EgArr, Jsc_vs_Eg] = calcJsc

figson = 1;

%load AM1.5
AM15_data = xlsread('AM15.xls');
AM15_data = AM15_data';
%AM15_data = [1,2,3,4; 1,2,3,4]

if figson == 1
    
figure(600)
plot(AM15_data(1, :), AM15_data(2, :));
xlabel('Wavelength [nm]');
ylabel('Power density [mWcm-2nm-1]');

end

h = 4.135667662e-15;        % eVs-1
c = 29979245800;            % cms-1
Vth = 0.0257;               % eV
q = 1.60217662e-19;         % C

t = 400e-7;
E = 1.6:0.01:100;

% Convert AM1.5 to correct units
E_ph = 1239.8./AM15_data(1, :);            % Photon energy
AM15 = AM15_data(2, :);     %(1e-3*(1/q)*

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