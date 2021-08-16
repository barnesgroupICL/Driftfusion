% An example script for comparing the current-voltage curves for two different
% architects of perovskite solar cell
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
tic
initialise_df

% Create parameters objects for Spiro/MAPI/TiO2 and PEDOT/MAPI/PCBM devices
par_tio2 = pc('Input_files/spiro_mapi_tio2.csv');
par_pcbm = pc('Input_files/pedotpss_mapi_pcbm.csv');

% Find equilibrium solutions
soleq_tio2 = equilibrate(par_tio2);
soleq_pcbm = equilibrate(par_pcbm);

%% Perform current-voltage scans
sol_CV_100mVs_tio2_dark = doCV(soleq_tio2.ion, 0, 0, 1.0, 0, 100e-3, 1, 241);
sol_CV_100mVs_pcbm_dark = doCV(soleq_pcbm.ion, 0, 0, 1.0, 0, 100e-3, 1, 241);
sol_CV_100mVs_tio2_1sun = doCV(soleq_tio2.ion, 1, 0, 1.2, 0, 100e-3, 1, 241);
sol_CV_100mVs_pcbm_1sun = doCV(soleq_pcbm.ion, 1, 0, 1.2, 0, 100e-3, 1, 241);

%% plot the current voltage curves
dfplot.JtotVapp(sol_CV_100mVs_tio2_dark,0)
hold on
dfplot.JtotVapp(sol_CV_100mVs_pcbm_dark,0)
hold on
dfplot.JtotVapp(sol_CV_100mVs_tio2_1sun,0)
hold on
dfplot.JtotVapp(sol_CV_100mVs_pcbm_1sun,0)
hold off

legend('tio2-dk', 'pcbm-dk', 'tio2-1sun', 'pcbm-1sun')
xlim([0, 1.2])
ylim([-30e-3, 10e-3])
hold off

%% plot the energy level diagram and carrier densities for the tio2 device at
% 1 V (t= 20s) during the illuminated forward scan
dfplot.ELxnpxacx(sol_CV_100mVs_tio2_1sun, 20)

% plot the currents as a function of position in the PCBM device at 0.5 V
% (t = 10 s) during the illuminated forward scan
dfplot.Jx(sol_CV_100mVs_tio2_1sun, 10)

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder sizes
% save('/Users/Username/Data/temp.mat')
toc
