%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Windows
% par_sio2 = pc('.\Input_files\pog2');
% Filepath Mac
par_sio2 = pc('Input_files/pog2.csv');
%% Equilibrium solutions
soleq_sio2 = equilibrate(par_sio2);

%% Plot equilibrium energy level diagram
dfplot.ELnpx(soleq_sio2.ion)

%% Current-voltage scan
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
% JVsol = doJV(soleq_sio2.ion, 100e-3, 201, 1, 0, 0, 1, 1);

% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV = doCV(soleq_sio2.ion, 1, 0, 1, 0, 100e-3, 1, 201);

%% Plot JV scan
dfplot.JtotVapp(sol_CV, 1);
set(gca,'YScale','log')

ylim([-30e-3,10e-3])
xlim([-0.2, 1.2])
