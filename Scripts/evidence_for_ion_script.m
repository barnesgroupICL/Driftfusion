%% A script to reproduce the results from Calado et al. 2016
% The results obtained here will be approximately the same as those
% presented in Calado 2016. Note that there may be some differences in the
% exact solutions owing to the difference in the methods used. Previously
% the Current-Voltage scan was reproduced using a discretised step-wise
% method (one solution per point) whereas here a single solution is used.
% Similarly the open circuit voltage transients were obtained using a
% symmetric model method previously while here a series resistance is
% applied to the boundary conditions.
% The recombination coefficients have also been adjusted as the finer mesh
% results in a larger SRH recombination rate close to the interfaces
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
initialise_df

%% Load parameters
% set OM =0 in pc before doing this for uniform generation in the active
% layer
% BC = Bottom Cathode i.e. the device with SRH recombination in the
% contacts (interfacial recombination)
% TC = Top Cathode i.e. device without SRH recombination
par_bc = pc('Input_files/Evidence_for_ion_migration_SRH.csv');
par_tc = pc('Input_files/Evidence_for_ion_migration_noSRH.csv');

% Get equilibrium solutions
soleq_bc = equilibrate(par_bc);
soleq_tc = equilibrate(par_tc);

%% Current-voltage scans

%% Bottom cathode device
% Scan to -1.0 V - dark only to get starting condition
% for the JV: doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
% Slow scan with ion mobility set to 1 cm2V-1s-1 (by setting mobseti =
% 1e10)- this ensure ions are in in equilibrium throughout scan

% Scan from -1.0 V to 1.2 V dark and light - ions must be accelerated by
%sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_bc = doCV(soleq_bc.ion, 1, 0, -1, 1.2, 40e-3, 1, 221);

%% Top cathode device
% Scan to -1.0 V with ion mobility off- dark only to get starting condition
% for the JV
sol_CV_tc = doCV(soleq_tc.ion, 1, 0, -1, 1.2, 40e-3, 1, 221);

%% Open circuit voltage tranisents
% Using Rs = 1 MOhm rather than the mirrored cell approach.
% lightonRs(sol_ini, Int, stab_time, mobseti, Rs, pnts)
Voc_transient_bc = lightonRs(soleq_bc.ion, 1, 100, 1, 1e6, 200);
Voc_transient_tc = lightonRs(soleq_tc.ion, 1, 100, 1, 1e6, 200);

%% Plot the JV's
dfplot.JtotVapp(sol_CV_bc, 0);
hold on
dfplot.JtotVapp(sol_CV_tc, 0);
xlim([-0.2,1])
ylim([-20e-3,5e-3])
hold off

%% Plot the Voc transients
dfplot.Voct(Voc_transient_bc);
hold on
dfplot.Voct(Voc_transient_tc);
hold off

% Plot the EL diagrams at 1 ms and 1 sec
dfplot.ELxnpxacx(Voc_transient_bc, [1e-3, 1])
% Uncomment to plot EL and carrier density diagrams for Top Cathode device
%dfplot.ELxnpxacx(Voc_transient.tc, [1e-3, 1])

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder sizes.
% save('/Users/Username/Data/temp.mat')
