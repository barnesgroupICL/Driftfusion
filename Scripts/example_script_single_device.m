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
par_tio2 = pc('Input_files/spiro_mapi_tio2.csv');

%% Find the equilibrium solutions
soleq_tio2 = equilibrate(par_tio2);

%% Perform dark and light current-voltage scan at 50 mVs-1 from 0 V to 1.2 V
% Input arguments: 
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JVsol_tio2 = doJV(soleq_tio2.ion, 50e-3, 100, 1, 1, 0, 1.2, 3);

%% plot the current voltage curve
dfplot.JV(JVsol_tio2,3)
ylim([-30e-3, 10e-3])

%% plot the energy level diagram and carrier densities for the tio2 device at
% 1 V (t= 20s) during the illuminated forward scan
dfplot.ELxnpxacx(JVsol_tio2.ill.f, 20)

%% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder sizes. This filepath may need to be altered for Windows-
% based systems
save('~/MATLAB_Data/temp.mat')
