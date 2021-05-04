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
par_3l = pc('Input_files/3_layer_test.csv');

%% Find the equilibrium solutions
soleq_3l = equilibrate(par_3l);

%% Perform dark and light current-voltage scan at 50 mVs-1 from 0 V to 1.2 V
% Input arguments: 
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
Vapp_func = 'sin';
DC_offset = 0;
Delta_AC = 0.02;
frequency = 1e3;
phase = 0;
Vapp_coeff = [DC_offset, Delta_AC, frequency, phase];

sol_sin = VappFunction(soleq_3l.ion, Vapp_func, Vapp_coeff, 1/frequency, 400, 0);

%% Plots
t = sol_sin.t;
Vapp = dfana.calcVapp(sol_sin);
J = dfana.calcJ(sol_sin);

figure(500)
plot(t, Vapp/max(Vapp), t, J.tot(:,1)/max(J.tot(:,1)))
xlabel('Time (s)')
ylabel('Vapp, J (norm.)')
%% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder sizes. This filepath may need to be altered for Windows-
% based systems
% save('~/MATLAB_Data/temp.mat')
