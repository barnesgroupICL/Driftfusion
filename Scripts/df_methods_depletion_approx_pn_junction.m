initialise_df

%% Create parameters objects
par_pn_base = pc('Input_files/pn_junction_DA_compare.csv');

soleq_base = equilibrate(par_pn_base, 1);

%% Base parameters current-voltage
% JV = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_base = doJV(soleq_base.el, 1e-9, 241, 1, 0, 0, 0.8, 3, 1);

%% Plot JV
dfplot.JV(JV_base,3);
hold on
ylim([-60e-3, 10e-3])

%% Get depletion approximation solution
[JV_DA.base, r0_DA.base, Voc_DA.base, g0_DA.base] = calcR0_pn_junction(par_pn_base);

%% tau = 1e-7 s
par_pn_tau1em7 = par_pn_base;
par_pn_tau1em7.taun(1) = 1e-7;
par_pn_tau1em7.taup(end) = 1e-7;
% Refresh the device to rebuild device structures
par_pn_tau1em7 = refresh_device(par_pn_tau1em7);

soleq_tau1em7 = equilibrate(par_pn_tau1em7, 1);

%% JV
JV_tau1em7 = doJV(soleq_tau1em7.el, 1e-9, 241, 1, 0, 0, 0.8, 3, 1);

%% Plots
dfplot.JV(JV_tau1em7,3);
hold on
ylim([-60e-3, 10e-3])

% Get depletion approximation solution
[JV_DA.tau1em7, r0_DA.tau1em7, Voc_DA.tau1em7, g0_DA.tau1em7] =  calcR0_pn_junction(par_pn_tau1em7);

%% tau = 1e-8 s
par_pn_tau1em8 = par_pn_base;
par_pn_tau1em8.taun(1) = 1e-8;
par_pn_tau1em8.taup(end) = 1e-8;
par_pn_tau1em8 = refresh_device(par_pn_tau1em8);

soleq_tau1em8 = equilibrate(par_pn_tau1em8, 1);

%% JV
JV_tau1em8 = doJV(soleq_tau1em8.el, 1e-9, 241, 1, 0, 0, 0.8, 3, 1);

%% Plots
dfplot.JV(JV_tau1em8,3);
hold on
ylim([-60e-3, 10e-3])

% Get depletion approximation solution
[JV_DA.tau1em8, r0_DA.tau1em8, Voc_DA.tau1em8, g0_DA.tau1em8] = calcR0_pn_junction(par_pn_tau1em8);

legend('DF, tau = 1e-6, dark', 'DF, tau = 1e-6, 1 sun', 'DA, tau = 1e-6, dark', 'DA, tau = 1e-6, 1 sun',...
    'DF, tau = 1e-7, dark', 'DF, tau = 1e-7, 1 sun', 'DA, tau = 1e-7, dark', 'DA, tau = 1e-7, 1 sun',...
    'DF, tau = 1e-8, dark', 'DF, tau = 1e-8, 1 sun', 'DA, tau = 1e-8, dark', 'DA, tau = 1e-8, 1 sun');