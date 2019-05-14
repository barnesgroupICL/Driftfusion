%% A script to reproduce the results from Calado et al. 2016

%% Load parameters
% set OM =0 in pc before doing this for uniform generation in the active
% layer
%par.tc = pc('input_files/Top Cathode.csv');
%par.bc = pc('input_files/Bottom Cathode.csv');

%% Current-voltage scans

%% Bottom cathode device
% Get equilibrium solution
%soleq.bc = equilibrate(par.bc);
% Scan to -1.0 V - dark only to get starting condition
% for the JV
% JV = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
%JV.bc_0V_to_m1V = doJV(soleq.bc.ion, 1e-6, 40, 1, 0, 0, -1, 1);
% Scan from -1.0 V to 1.4 V dark and light
JV.bc_m1V_to_0V = doJV(JV.bc_0V_to_m1V.dk.f, 40e-1, 200, 1, 1, -1, 0, 2);
JV.bc_0V_to_1p2V = doJV(JV.bc_m1V_to_0V.ill.f, 40e-1, 200, 1, 1, 0, 1.2, 2);

%% Top cathode device
% Get equilibrium solution
soleq.tc = equilibrate(par.tc);
% Scan to -1.0 V with ion mobility off- dark only to get starting condition
% for the JV
JV.tc_0V_to_m1V = doJV(soleq.tc.ion, 1e-6, 40, 1, 0, 0, -1, 1);
% Scan from -1.0 V to 1.4 V dark and light
JV.tc_m1V_to_0V = doJV(JV.tc_0V_to_m1V.dk.f, 40e-1, 200, 1, 1, -1, 0, 2);
JV.tc_0V_to_1p2V = doJV(JV.tc_m1V_to_0V.ill.f, 40e-1, 200, 1, 1, 0, 1.2, 2);

