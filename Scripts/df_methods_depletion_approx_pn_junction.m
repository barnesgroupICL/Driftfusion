initialise_df

%% Import pn junction properties from external .csv file
par_pn = pc('Input_files/pn_junction_DA_compare.csv');

%% Obtain equilibrium solution
soleq_pn = equilibrate(par_pn);

%% Current voltage scans
% JV = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_pn = doJV(soleq_pn.el, 1, 241, 1, 0, 0, 1.1, 3);

%% Plot equilibrium
dfplot.ELnpx(soleq_pn.el);