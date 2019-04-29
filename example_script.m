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

% plot the outputs
dfplot.JV(JV.tio2.JV50mVs,3)
figure(4)
hold on
dfplot.JV(JV.pcbm.JV50mVs,3)
legend('tio2-dk-f', 'tio2-dk-r', 'tio2-1sun-f', 'tio2-1sun-r',...
        'pcbm-dk-f', 'pcbm-dk-r', 'pcbm-1sun-f', 'pcbm-1sun-r')
hold off