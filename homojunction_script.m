%% Ionic screening in perovskite p-n homojunctions
% A script to reproduce the results from Calado, P. & Barnes, P. R. F. 
% Ionic screening in perovskite p?n homojunctions. Nat. Energy 13?15 (2021). 
% doi:10.1038/s41560-021-00838-1
% Script author: P Calado, 2021
% Email: p.calado13@imperial.ac.uk

%% Add folders and set plotting defaults
initialise_df

%% Intrinsic, Nion= 1e19 cm-3
%% Intrinsic
% Create a parameters object for intrinsic absorber device
par_intrinsic = pc('input_files/intrinsic.csv');
% Find the equilibrium solutions with electronic carriers only, SOLEQ.EL
% and with electronic and ionic carriers SOLEQ.ION
soleq_intrinsic = equilibrate(par_intrinsic);

% Run current density voltage scans (J-V) at 1 mVs-1 from 0 to 1.2 V (the
% voltage is negative in the argument due to the direction in which the
% device is built in the simulation). Where _ION is used, this indicates
% that Ncat = 1e19 cm-3 is used.
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_intrinsic_el = doJV(soleq_intrinsic.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV_intrinsic_ion = doJV(soleq_intrinsic.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);

% Find approximate open circuit, illuminated steady-state solution.
% Inputting a negative time uses the modulus as the initial guess but loops
% with increasing time to find steady-state solution. Here we use a series
% resistance of 10^6 Ohm cm2 to approximate open circuit.
% sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
OC_intrinsic_el = lightonRs(soleq_intrinsic.el, 1, 1, 0, 1e6, 200);
OC_intrinsic_ion = lightonRs(soleq_intrinsic.ion, 1, -100, 1, 1e6, 200);

%% Intrinsic, Nion= 1e17 cm-3
par_intrinsic_Ncat1e17 = pc('input_files/intrinsic.csv');
% Here we overwrite the ionic densities for each layer. Although ions are
% defined in every layer, in the transport layers and interfaces,
% the mobilities are set to zero and the charge is balanced by a static counter ion density 
% such that the net space charge charge density due to ionic variables is zero.
par_intrinsic_Ncat1e17.Ncat = [1e17,1e17,1e17,1e17,1e17,1e17,1e17];
par_intrinsic_Ncat1e17.Nani = [1e17,1e17,1e17,1e17,1e17,1e17,1e17];
par_intrinsic_Ncat1e17 = refresh_device(par_intrinsic_Ncat1e17);

soleq_intrinsic_Ncat1e17 = equilibrate(par_intrinsic_Ncat1e17);

% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_intrinsic_Ncat1e17_ion = doJV(soleq_intrinsic_Ncat1e17.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
% sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
OC_intrinsic_Ncat1e17_ion = lightonRs(soleq_intrinsic_Ncat1e17.ion, 1, -100, 1, 1e6, 200);

%% Intrinsic, Nion= 1e18 cm-3
par_intrinsic_Ncat1e18 = pc('input_files/intrinsic.csv');
par_intrinsic_Ncat1e18.Ncat = [1e18,1e18,1e18,1e18,1e18,1e18,1e18];
par_intrinsic_Ncat1e18.Nani = [1e18,1e18,1e18,1e18,1e18,1e18,1e18];
par_intrinsic_Ncat1e18 = refresh_device(par_intrinsic_Ncat1e18);

soleq_intrinsic_Ncat1e18 = equilibrate(par_intrinsic_Ncat1e18);

% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_intrinsic_Ncat1e18_ion = doJV(soleq_intrinsic_Ncat1e18.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
% sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
OC_intrinsic_Ncat1e18_ion = lightonRs(soleq_intrinsic_Ncat1e18.ion, 1, -100, 1, 1e6, 200);

%% n-type
par_ntype = pc('input_files/n_type.csv');
soleq_ntype = equilibrate(par_ntype);

% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_ntype_el = doJV(soleq_ntype.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV_ntype_ion = doJV(soleq_ntype.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
% sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
OC_ntype_el = lightonRs(soleq_ntype.el, 1, 1, 0, 1e6, 200);
OC_ntype_ion = lightonRs(soleq_ntype.ion, 1, -100, 0, 1e6, 200);

%% Homojunction
par_hom = pc('input_files/homojunction.csv');
soleq_hom = equilibrate(par_hom);

% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_hom_el = doJV(soleq_hom.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV_hom_ion = doJV(soleq_hom.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
% sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
OC_hom_el = lightonRs(soleq_hom.el, 1, 1, 0, 1e6, 200);
OC_hom_ion = lightonRs(soleq_hom.ion, 1, -100, 1, 1e6, 200);

%% Homojunction, Nion= 1e17 cm-3
par_hom_Ncat1e17 = pc('input_files/homojunction.csv');
par_hom_Ncat1e17.Ncat = [1e17,1e17,1e17,1e17,1e17,1e17,1e17];
par_hom_Ncat1e17.Nani = [1e17,1e17,1e17,1e17,1e17,1e17,1e17];
par_hom_Ncat1e17 = refresh_device(par_hom_Ncat1e17);

soleq_hom_Ncat1e17 = equilibrate(par_hom_Ncat1e17);
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_hom_Ncat1e17_ion = doJV(soleq_hom_Ncat1e17.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
% sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
OC_hom_Ncat1e17_ion = lightonRs(soleq_hom_Ncat1e17.ion, 1, -100, 1, 1e6, 200);

%% Homojunction, Nion= 1e18 cm-3
par_hom_Ncat1e18 = pc('input_files/homojunction.csv');
par_hom_Ncat1e18.Ncat = [1e18,1e18,1e18,1e18,1e18,1e18,1e18];
par_hom_Ncat1e18.Nani = [1e18,1e18,1e18,1e18,1e18,1e18,1e18];
par_hom_Ncat1e18 = refresh_device(par_hom_Ncat1e18);

soleq_hom_Ncat1e18 = equilibrate(par_hom_Ncat1e18);
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_hom_Ncat1e18_ion = doJV(soleq_hom_Ncat1e18.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
% sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
OC_hom_Ncat1e18_ion = lightonRs(soleq_hom_Ncat1e18.ion, 1, -100, 1, 1e6, 200);

%% J turn ons, Nion = 1e19 cm-3
Jmpp_intrinsic_el = jumptoV(soleq_intrinsic.el, -0.85, 10, 0, 1, 0, 0);
Jmpp_intrinsic_ion = jumptoV(soleq_intrinsic.ion, -0.92, 10, 1, 1, 0, 0);
Jmpp_hom_el = jumptoV(soleq_hom.el, -0.93, 10, 0, 1, 0, 0);
Jmpp_hom_ion = jumptoV(soleq_hom.ion, -0.95, 10, 1, 1, 0, 0);

%% J-V 1d Figure from Calado & Barnes (forward scan indicated by -F, reverse scan indicated by -R)
%% Mobile electronic and ionic carriers
dfplot.JV(JV_intrinsic_el, 2)
figure(4)
hold on
dfplot.JV(JV_intrinsic_ion, 2)
figure(4)
hold on
dfplot.JV(JV_hom_el, 2)
figure(4)
hold on
dfplot.JV(JV_hom_ion, 2)
figure(4)
hold on
dfplot.JV(JV_ntype_el, 2)
figure(4)
hold off
legend('intrinsic -F', 'intrinsic -R', 'intrinsic, Ncat = 1e19 cm-3 -F', 'intrinsic, Ncat = 1e19 cm-3 -R',...
    'homojunciton - F', 'homojunciton - R', 'homojunciton, Ncat = 1e19 cm-3 - F', 'homojunciton, Ncat = 1e19 cm-3 - R',...
    'n-type -F', 'n-type -R')
xlim([0.7, 1.1])
ylim([-5e-3, 25e-3])

%% Open circuit Electric field
dfplot.Fx(OC_hom_el)
hold on
dfplot.Fx(OC_hom_Ncat1e17_ion)
hold on
dfplot.Fx(OC_hom_Ncat1e18_ion)
hold on
dfplot.Fx(OC_hom_ion)           % Nion = 1e19 cm-3
legend('HTL','int1','n-type','int2','p-type','int3','ETL',...
    'Homojunction, N_{cat} = 0 cm{-3}', 'Homojunction, N_{cat} = 10^{17} cm{-3}',...
    'Homojunction, N_{cat} = 10^{18} cm{-3}', 'Homojunction, N_{cat} = 10^{19} cm{-3}')
hold off
ylim([-0.5e5,2e5])

%% Steady state maximum power current density transients
dfplot.Jtott(Jmpp_intrinsic_el, Jmpp_intrinsic_el.x(end));
hold on
dfplot.Jtott(Jmpp_intrinsic_ion, Jmpp_intrinsic_ion.x(end));
hold on
dfplot.Jtott(Jmpp_hom_el, Jmpp_hom_el.x(end));
hold on
dfplot.Jtott(Jmpp_hom_ion, Jmpp_hom_ion.x(end));
hold off
ylim([0.0215, 0.022])

%% Short circuit energy level diagrams and carriers for the intrinsic device
%% Electronic carriers only
dfplot.ELxnpx(soleq_intrinsic.el)
subplot(2,1,1)
hold on
subplot(2,1,2);
hold on

dfplot.ELxnpx(soleq_ntype.el)
subplot(2,1,1)
ylim([-8, -1])
legend('', '', '', '', '',... % Layer blocks
    'intrinsic, \it{E}_{fn}', 'intrinsic, \it{E}_{fp}', 'intrinsic, \it{E}_{CB}', 'intrinsic, \it{E}_{VB}',...
    'n-type, \it{E}_{fn}', 'n-type, \it{E}_{fp}', 'n-type, \it{E}_{CB}', 'n-type, \it{E}_{VB}')

subplot(2,1,2);
ylim([1e-4,1e20])
legend('', '', '', '', '',... % Layer blocks
    'intrinsic, \it{n}', 'intrinsic, \it{p}',...
    'n-type, \it{n}', 'n-type, \it{p}')


%% Electronic and ionic carriers
dfplot.ELxnpxacx(soleq_intrinsic_Ncat1e17.ion)
subplot(3,1,1)
subplot(3,1,2);
hold on
subplot(3,1,3);
hold on

dfplot.ELxnpxacx(soleq_intrinsic_Ncat1e18.ion)
subplot(3,1,1)
subplot(3,1,2);
hold on
subplot(3,1,3);
hold on

dfplot.ELxnpxacx(soleq_intrinsic.ion)
subplot(3,1,1)
ylim([-8, -1])
subplot(3,1,2);
legend('', '', '', '', '', '', '',... % Layer blocks
    '\it{N}_{cat} = 10^{17} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{17} cm^{-3}, \it{p}',...
    '\it{N}_{cat} = 10^{18} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{18} cm^{-3}, \it{p}',...
    '\it{N}_{cat} = 10^{19} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{19} cm^{-3}, \it{p}')
hold off
ylim([1e10,1e20])
subplot(3,1,3);
legend('', '', '', '', '', '', '',... % Layer blocks
    '\it{N}_{cat} = 10^{17} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{17} cm^{-3}, \it{a}',...
    '\it{N}_{cat} = 10^{18} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{18} cm^{-3}, \it{a}',...
    '\it{N}_{cat} = 10^{19} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{19} cm^{-3}, \it{a}')
hold off
ylim([0, 4e19])

%% Open circuit energy level diagrams and carriers for the intrinsic device
%% Electronic carriers only
dfplot.ELxnpx(OC_intrinsic_el)
subplot(2,1,1)
hold on
subplot(2,1,2);
hold on

dfplot.ELxnpx(OC_ntype_el)
subplot(2,1,1)
ylim([-8, -2])
legend('', '', '', '', '',... % Layer blocks
    'intrinsic, \it{E}_{fn}', 'intrinsic, \it{E}_{fp}', 'intrinsic, \it{E}_{CB}', 'intrinsic, \it{E}_{VB}',...
    'n-type, \it{E}_{fn}', 'n-type, \it{E}_{fp}', 'n-type, \it{E}_{CB}', 'n-type, \it{E}_{VB}')

subplot(2,1,2);
ylim([1e14,1e20])
legend('', '', '', '', '',... % Layer blocks
    'intrinsic, \it{n}', 'intrinsic, \it{p}',...
    'n-type, \it{n}', 'n-type, \it{p}')


%% Electronic and ionic carriers
dfplot.ELxnpxacx(OC_intrinsic_Ncat1e17_ion)
subplot(3,1,1)
subplot(3,1,2);
hold on
subplot(3,1,3);
hold on

dfplot.ELxnpxacx(OC_intrinsic_Ncat1e18_ion)
subplot(3,1,1)
subplot(3,1,2);
hold on
subplot(3,1,3);
hold on

dfplot.ELxnpxacx(OC_intrinsic_ion)
subplot(3,1,1)
ylim([-8, -1])
subplot(3,1,2);
legend('', '', '', '', '', '', '',... % Layer blocks
    '\it{N}_{cat} = 10^{17} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{17} cm^{-3}, \it{p}',...
    '\it{N}_{cat} = 10^{18} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{18} cm^{-3}, \it{p}',...
    '\it{N}_{cat} = 10^{19} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{19} cm^{-3}, \it{p}')
hold off
ylim([1e14,1e20])
subplot(3,1,3);
legend('', '', '', '', '', '', '',... % Layer blocks
    '\it{N}_{cat} = 10^{17} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{17} cm^{-3}, \it{a}',...
    '\it{N}_{cat} = 10^{18} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{18} cm^{-3}, \it{a}',...
    '\it{N}_{cat} = 10^{19} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{19} cm^{-3}, \it{a}')
hold off
ylim([0, 1.2e19])

%% Short circuit energy level diagrams and carriers for homojunction
%% Electronic carriers only
dfplot.ELxnpxacx(soleq_hom.el)
subplot(3,1,1)
ylim([-8, -1])
subplot(3,1,2);
ylim([1e-4,1e20])
subplot(3,1,3);
ylim([0, 4e19])
%% Electronic and ionic carriers
dfplot.ELxnpxacx(soleq_hom_Ncat1e17.ion)
subplot(3,1,1)
subplot(3,1,2);
hold on
subplot(3,1,3);
hold on

dfplot.ELxnpxacx(soleq_hom_Ncat1e18.ion)
subplot(3,1,1)
subplot(3,1,2);
hold on
subplot(3,1,3);
hold on

dfplot.ELxnpxacx(soleq_hom.ion)
subplot(3,1,1)
ylim([-8, -1])
subplot(3,1,2);
legend('', '', '', '', '', '', '',... % Layer blocks
    '\it{N}_{cat} = 10^{17} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{17} cm^{-3}, \it{p}',...
    '\it{N}_{cat} = 10^{18} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{18} cm^{-3}, \it{p}',...
    '\it{N}_{cat} = 10^{19} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{19} cm^{-3}, \it{p}')
hold off
ylim([1e10,1e20])
subplot(3,1,3);
legend('', '', '', '', '', '', '',... % Layer blocks
    '\it{N}_{cat} = 10^{17} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{17} cm^{-3}, \it{a}',...
    '\it{N}_{cat} = 10^{18} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{18} cm^{-3}, \it{a}',...
    '\it{N}_{cat} = 10^{19} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{19} cm^{-3}, \it{a}')
hold off
ylim([0, 4e19])

%% Open circuit energy level diagrams and carriers for homojunction
%% Electronic carriers only
dfplot.ELxnpxacx(OC_hom_el)
subplot(3,1,1)
ylim([-8, -2])
subplot(3,1,2)
ylim([1e14,1e20])
subplot(3,1,3);
ylim([0, 1.2e19])

%% Electronic and ionic carriers
dfplot.ELxnpxacx(OC_hom_Ncat1e17_ion)
subplot(3,1,1)
subplot(3,1,2);
hold on
subplot(3,1,3);
hold on

dfplot.ELxnpxacx(OC_hom_Ncat1e18_ion)
subplot(3,1,1)
subplot(3,1,2);
hold on
subplot(3,1,3);
hold on

dfplot.ELxnpxacx(OC_hom_ion)
subplot(3,1,1)
ylim([-8, -2])
subplot(3,1,2);
legend('', '', '', '', '', '', '',... % Layer blocks
    '\it{N}_{cat} = 10^{17} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{17} cm^{-3}, \it{p}',...
    '\it{N}_{cat} = 10^{18} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{18} cm^{-3}, \it{p}',...
    '\it{N}_{cat} = 10^{19} cm^{-3}, \it{n}', '\it{N}_{cat} = 10^{19} cm^{-3}, \it{p}')
hold off
ylim([1e14,1e20])
subplot(3,1,3);
legend('', '', '', '', '', '', '',... % Layer blocks
    '\it{N}_{cat} = 10^{17} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{17} cm^{-3}, \it{a}',...
    '\it{N}_{cat} = 10^{18} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{18} cm^{-3}, \it{a}',...
    '\it{N}_{cat} = 10^{19} cm^{-3}, \it{c}', '\it{N}_{cat} = 10^{19} cm^{-3}, \it{a}')
hold off
ylim([0, 1.2e19])