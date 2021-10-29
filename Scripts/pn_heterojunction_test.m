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
par_pn_hetero = pc('Input_files/pn_heterojunction.csv');

%% Find the equilibrium solutions
soleq_pn_hetero = equilibrate(par_pn_hetero);

%% Perform dark and light current-voltage scan at 50 mVs-1 from 0 V to 1.2 V
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_100mVs_pn_hetero = doCV(soleq_pn_hetero.el, 0, 0, 0.8, -0.2, 100e-3, 1, 281);

%% plot the energy level diagram and carrier densities for the device at
% 1 V (t= 10s) during the illuminated forward scan
dfplot.ELxnpxacx(sol_CV_100mVs_pn_hetero, 10)

%% plot J-V
dfplot.JVapp(sol_CV_100mVs_pn_hetero, par_pn_hetero.d_midactive)
%ylim([-20e-3, 20e-3])

%% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder sizes. This filepath may need to be altered for Windows-
% based systems
% save('~/MATLAB_Data/temp.mat')
