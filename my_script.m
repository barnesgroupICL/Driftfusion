%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
par_sio2 = pc('C:\Users\Manan Mehta\Documents\GitHub\Driftfusion\Input_files\pog2');

%% Equilibrium solutions
soleq_sio2 = equilibrate(par_sio2);
 %% Equilibrium solutions
JVsol = doJV(soleq_sio2.ion, 100e-3, 201, 1, 0, 0, 1, 1);

% plot JV scan
dfplot.JV(JVsol, 1);
set(gca,'YScale','log')

ylim([-30e-3,10e-3])
xlim([-0.2, 1.2])


