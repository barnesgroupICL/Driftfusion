initialise_df

% Create parameters objects for Spiro/MAPI/TiO2 and PEDOT/MAPI/PCBM devices
par.tio2 = pc('input_files/tio2.csv');
par.pcbm = pc('input_files/pcbm.csv');

% Find equilibrium solutions
soleq.tio2 = equilibrate(par.tio2);
soleq.pcbm = equilibrate(par.pcbm);

% Perform current-voltage scans to 1.2 V
JV.tio2.JV50mVs = doJV(soleq.tio2.i_sr, 50e-3, 100, 1, 1, 0, 1.2, 3);
JV.pcbm.JV50mVs = doJV(soleq.pcbm.i_sr, 50e-3, 100, 1, 1, 0, 1.2, 3);

% plot the current voltage curve
dfplot.JV(JV.tio2.JV50mVs,3)

figure(4)
hold on

dfplot.JV(JV.pcbm.JV50mVs,3)
legend('tio2-dk-f', 'tio2-dk-r', 'tio2-1sun-f', 'tio2-1sun-r',...
        'pcbm-dk-f', 'pcbm-dk-r', 'pcbm-1sun-f', 'pcbm-1sun-r')
hold off

% plot the energy level diagram and charge densities for the tio2 device at 
% 1 V (t= 20s) during the illuminated forward scan
dfplot.ELx(JV.tio2.JV50mVs.ill.f, 20)

% plot the currents as a function of position in the PCBM device at 0.5 V
% (t = 10 s) during the illuminated forward scan
dfplot.Jx(JV.pcbm.JV50mVs.ill.f, 10)

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder size regardless of whether it contains data or not.
% save('/Users/Username/Data/temp.mat')
