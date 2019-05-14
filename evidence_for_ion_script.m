%% A script to reproduce the results from Calado et al. 2016
% The results obtained here will be approximately the same as those
% presented in Calado 2016. Note that there may be some differences in the
% exact solutions owing to the difference in the methods used. Previously
% the Current-Voltage scan was reproduced using a discretised step wise
% method (one solution per point) whereas here a single solution is used.
% Similarly the open circuit voltage transients were obtained using a
% symmetric model method previously while here a series resistance is
% applied to the boundary conditions.

%% Load parameters
% set OM =0 in pc before doing this for uniform generation in the active
% layer
par.tc = pc('input_files/Top Cathode.csv');
par.bc = pc('input_files/Bottom Cathode.csv');

% Get equilibrium solutions
soleq.bc = equilibrate(par.bc);
soleq.tc = equilibrate(par.tc);

%% Current-voltage scans

%% Bottom cathode device
% Scan to -1.0 V - dark only to get starting condition
% for the JV: doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
% Slow scan with ion mobility set to 1 cm2V-1s-1 (by setting mobseti =
% 1e10)- this ensure ions are in in equilibrium throughout scan
JV.bc_0V_to_m1V = doJV(soleq.bc.ion, 1e-6, 40, 1, 1e10, 0, -1, 1);
% Scan from -1.0 V to 1.2 V dark and light
JV.bc_m1V_to_1p2V = doJV(JV.bc_0V_to_m1V.dk.f, 40e-1, 200, 1, 1, -1, 1.2, 3);

% Set JV plot limits
figure(4)
xlim([-0.2,1])
ylim([-20e-3,5e-3])
hold on

%% Top cathode device
% Scan to -1.0 V with ion mobility off- dark only to get starting condition
% for the JV
JV.tc_0V_to_m1V = doJV(soleq.tc.ion, 1e-6, 40, 1, 1e10, 0, -1, 1);
% Scan from -1.0 V to 1.4 V dark and light
JV.tc_m1V_to_1p2V = doJV(JV.tc_0V_to_m1V.dk.f, 40e-1, 200, 1, 1, -1, 1.2, 3);

figure(4)
hold off

%% Open circuit voltage tranisents
% lighton_Rs(sol_ini, Int, stab_time, mobseti, Rs, pnts)
Voc_transient.bc = lighton_Rs(soleq.bc.ion, 1, 1, 1, 1e6, 400);
Voc_transient.tc = lighton_Rs(soleq.tc.ion, 1, 1, 1, 1e6, 400);
