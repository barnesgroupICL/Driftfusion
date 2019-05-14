%% A script to reproduce the results from Calado et al. 2016
% The results obtained here will be approximately the same as those
% presented in Calado 2016. Note that there may be some differences in the
% exact solutions owing to the difference in the methods used. Previously
% the Current-Voltage scan was reproduced using a discretised step-wise
% method (one solution per point) whereas here a single solution is used.
% Similarly the open circuit voltage transients were obtained using a
% symmetric model method previously while here a series resistance is
% applied to the boundary conditions.
% Note also that the ions have been accelerated by two orders of magnitude
% as compared to the original publication mu_ion = 10^-10 cm2V-1s-1 to aid 
% convergence of the solutions

%% Load parameters
% set OM =0 in pc before doing this for uniform generation in the active
% layer
% BC = Bottom Cathode i.e. the device with SRH recombination in the
% contacts (interfacial recombination)
% TC = Top Cathode i.e. device without SRH recombination
par.tc = pc('input_files/Evidence_for_ion_migration_noSRH.csv');
par.bc = pc('input_files/Evidence_for_ion_migration_SRH.csv');

% Get equilibrium solutions
soleq.bc = equilibrate(par.bc);
soleq.tc = equilibrate(par.tc);

%% Current-voltage scans

%% Bottom cathode device
% Scan to -1.0 V - dark only to get starting condition
% for the JV: doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
% Slow scan with ion mobility set to 1 cm2V-1s-1 (by setting mobseti =
% 1e10)- this ensure ions are in in equilibrium throughout scan
JV.bc_0V_to_m1V = doJV(soleq.bc.ion, 1e-6, 40, 0, 1e10, 0, -1, 1);
% Scan from -1.0 V to 1.2 V dark and light - ions must be accelerated by
% x10 for convergence
JV.bc_m1V_to_1p2V = doJV(JV.bc_0V_to_m1V.dk.f, 40, 200, 1, 10, -1, 1.2, 3);

%% Top cathode device
% Scan to -1.0 V with ion mobility off- dark only to get starting condition
% for the JV
JV.tc_0V_to_m1V = doJV(soleq.tc.ion, 1e-6, 40, 0, 1e10, 0, -1, 1);
% Scan from -1.0 V to 1.4 V dark and light
JV.tc_m1V_to_1p2V = doJV(JV.tc_0V_to_m1V.dk.f, 40, 200, 1, 10, -1, 1.2, 3);

% Plot the JV's and set limits
dfplot.JV(JV.bc_m1V_to_1p2V,3);
hold on
dfplot.JV(JV.tc_m1V_to_1p2V,3);
xlim([-0.2,1])
ylim([-20e-3,5e-3])
hold off

%% Open circuit voltage tranisents
% Using Rs = 1 MOhm rather than the mirrored cell approach.
% lighton_Rs(sol_ini, Int, stab_time, mobseti, Rs, pnts)
Voc_transient.bc = lighton_Rs(soleq.bc.ion, 1, 1, 1, 1e6, 400);
Voc_transient.tc = lighton_Rs(soleq.tc.ion, 1, 1, 1, 1e6, 400);

% Plot the Voc transients
dfplot.Voct(Voc_transient.bc);
hold on
dfplot.Voct(Voc_transient.tc);
hold off

% Plot the EL diagrams at 1 ms and 1 sec
dfplot.ELx(Voc_transient.bc, [1e-3, 1])
% Uncomment to plot EL and carrier density diagrams for Top Cathode device
%dfplot.ELx(Voc_transient.tc, [1e-3, 1])