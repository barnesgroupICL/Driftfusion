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
% Initialise the system
initialise_df

%% Create a parameters object for Spiro/MAPI/TiO2 by including a filepath to the 
% appropriate .csv as the arugment to the parameters class PC
par_pedot_mapi_pcbm = pc('Input_files/pedotpss_mapi_pcbm.csv');

%% Find the equilibrium solutions
soleq_pedot_mapi_pcbm = equilibrate(par_pedot_mapi_pcbm);

%% Perform dark and light current-voltage scan at 50 mVs-1 from 0 V to 1.2 V
% Input arguments: 
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
sol_CV_100mVs_pedot_mapi_pcbm = doCV(soleq_pedot_mapi_pcbm.ion, 1, 0, 1.2, -0.2, 100e-3, 1, 281);
sol_CV_200mVs_pedot_mapi_pcbm = doCV(soleq_pedot_mapi_pcbm.ion, 1, 0, 1.2, -0.2, 200e-3, 1, 281);
sol_CV_400mVs_pedot_mapi_pcbm = doCV(soleq_pedot_mapi_pcbm.ion, 1, 0, 1.2, -0.2, 400e-3, 1, 281);

%% plot the current voltage curve
dfplot.JtotVapp(sol_CV_100mVs_pedot_mapi_pcbm, 0)
hold on
dfplot.JtotVapp(sol_CV_200mVs_pedot_mapi_pcbm, 0)
hold on
dfplot.JtotVapp(sol_CV_400mVs_pedot_mapi_pcbm, 0)
hold off
ylim([-30e-3,10e-3])
xlim([-0.2, 1.2])
legend('100 mVs-1', '200 mVs-1', '400 mVs-1') 

%% plot the energy level diagram and carrier densities for the tio2 device at
% 1 V (t= 10s) during the illuminated forward scan
dfplot.ELnpx(sol_CV_100mVs, 10)
%% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder sizes. This filepath may need to be altered for Windows-
% based systems
save('~/MATLAB_Data/temp.mat')
