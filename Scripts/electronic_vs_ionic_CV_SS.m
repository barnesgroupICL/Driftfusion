%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%% Start code
% Initialise the system
initialise_df

%% Create a parameters object for Spiro/MAPI/TiO2 by including a filepath to the 
% appropriate .csv as the arugment to the parameters class PC
par_tio2 = pc('./Input_files/spiro_mapi_tio2.csv');

%% Find the equilibrium solutions
soleq_tio2 = equilibrate(par_tio2);

%% Perform dark and light current-voltage scan at 50 mVs-1 from 0 V to 1.2 V
% Input arguments: 
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
% sol_CV_100mVs_dark = doCV(soleq_tio2.ion, 0, 0, 1.2, -0.2, 100e-3, 1, 281);
% sol_CV_200mVs_dark = doCV(soleq_tio2.ion, 0, 0, 1.2, -0.2, 200e-3, 1, 281);
% sol_CV_400mVs_dark = doCV(soleq_tio2.ion, 0, 0, 1.2, -0.2, 400e-3, 1, 281);
sol_CV_el_dk = doCV(soleq_tio2.el, 0, 0, 1.0, 0, 100e-3, 1, 201);
sol_CV_el = doCV(soleq_tio2.el, 1, 0, 1.0, 0, 100e-3, 1, 201);
%% With ions
soleq_tio2.ion.par.K_c = 1e3;
sol_CV_ion_dk = doCV(soleq_tio2.ion, 0, 0, 1.0, 0, 100e-3, 1, 201);
sol_CV_ion = doCV(soleq_tio2.ion, 1, 0, 1.0, 0, 100e-3, 1, 201);

%% plot the current voltage curve
Vbi_yplot = -50e-3:1e-3:50e-3;
Vbi_xplot = 0.6*ones(1, length(Vbi_yplot));

figure(1)
dfplot.JtotVapp(sol_CV_el, 0)
hold on
dfplot.JtotVapp(sol_CV_ion, 0)
hold on
plot(Vbi_xplot, Vbi_yplot, 'k--')
t0 = annotation('textbox', [0.6, 0.1, 0.1, 0.1], 'String', 'V_{bi}');
t0.FontSize = 18;
t0.EdgeColor = 'none';
hold off
xlabel('Applied voltage [V]')
ylabel('Current density [Acm^{-2}]')
ylim([-25e-3, 10e-3])
xlim([0, 1.0])
axes('Position',[.25 .50 .4 .4])
box on
dfplot.JtotVapp(sol_CV_el_dk, 0)
hold on
dfplot.JtotVapp(sol_CV_ion_dk, 0)
hold off
xlabel('Applied voltage [V]')
ylabel('Current density [Acm^{-2}]')
lgd = legend('N_{ion} = 0 cm^{-3}', 'N_{ion} = 10^{19} cm^{-3}');
lgd.EdgeColor = 'none'
set(gca, 'YScale', 'log')
set(gca,'FontSize',14)

%%
figure(2)
subplot(2, 1, 1)
dfplot.Vx(sol_CV_el_dk, [0, 2, 4, 6, 8])
t1 = annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String', 'N_{ion} = 0 cm^{-3}');
t1.FontSize = 14;
t1.EdgeColor = 'none';
t1a = annotation('textbox', [0.15, 0.52, 0.1, 0.1], 'String', 'HTL');
t1a.FontSize = 14;
t1a.EdgeColor = 'none';
t1b = annotation('textbox', [0.5, 0.52, 0.1, 0.1], 'String', 'Absorber');
t1b.FontSize = 14;
t1b.EdgeColor = 'none';
t1c = annotation('textbox', [0.8, 0.52, 0.1, 0.1], 'String', 'ETL');
t1c.FontSize = 14;
t1c.EdgeColor = 'none';
subplot(2, 1, 2)
dfplot.Vx(sol_CV_ion_dk, [0, 2, 4, 6, 8])
t2 = annotation('textbox', [0.15, 0.35, 0.1, 0.1], 'String', 'N_{ion} = 10^{19} cm^{-3}');
t2.FontSize = 14;
t2.EdgeColor = 'none';
t2a = annotation('textbox', [0.15, 0.05, 0.1, 0.1], 'String', 'HTL');
t2a.FontSize = 14;
t2a.EdgeColor = 'none';
t2b = annotation('textbox', [0.5, 0.05, 0.1, 0.1], 'String', 'Absorber');
t2b.FontSize = 14;
t2b.EdgeColor = 'none';
t2c = annotation('textbox', [0.8, 0.05, 0.1, 0.1], 'String', 'ETL');
t2c.FontSize = 14;
t2c.EdgeColor = 'none';
hold off
%% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder sizes. This filepath may need to be altered for Windows-
% based systems
% save('~/MATLAB_Data/temp.mat')
