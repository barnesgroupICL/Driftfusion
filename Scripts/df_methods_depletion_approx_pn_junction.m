initialise_df

%% Create parameters objects
% par = pc('Input_files/spiro_mapi_tio2_simple.csv');
par_pn_base = pc('Input_files/pn_junction_DA_compare.csv');

soleq_base = equilibrate_short(par_pn_base);

% JV = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_base = doJV(soleq_base, 1e-9, 241, 1, 0, 0, 0.8, 3);

dfplot.JV(JV_base,3);
hold on
ylim([-60e-3, 10e-3])

% Get depletion approximation solution
[JV_DA.base, r0_DA.base, Voc_DA.base, g0_DA.base] = calcR0_pn_junction(par_pn_base);

%% tau = 1e-5 s
par_pn_tau1em5 = par_pn_base;

par_pn_tau1em5.taun(1) = 1e-5;
par_pn_tau1em5.taup(end) = 1e-5;

par_pn_tau1em5 = refresh_device(par_pn_tau1em5);

soleq_tau1em5 = equilibrate_short(par_pn_tau1em5);

JV_tau1em5 = doJV_short(soleq_tau1em5, 1e-9, 241, 1, 0, 0, 0.8, 3);

dfplot.JV(JV_tau1em5,3);
hold on
ylim([-60e-3, 10e-3])

% Get depletion approximation solution
[JV_DA.tau1em5, r0_DA.tau1em5, Voc_DA.tau1em5, g0_DA.tau1em5] = calcR0_pn_junction(par_pn_tau1em5);

%% tau = 1e-7 s
par_pn_tau1em7 = par_pn_base;

par_pn_tau1em7.taun(1) = 1e-7;
par_pn_tau1em7.taup(end) = 1e-7;

par_pn_tau1em7 = refresh_device(par_pn_tau1em7);

soleq_tau1em7 = equilibrate_short(par_pn_tau1em7);

JV_tau1em7 = doJV_short(soleq_tau1em7, 1e-9, 241, 1, 0, 0, 0.8, 3);

dfplot.JV(JV_tau1em7,3);
hold on
ylim([-60e-3, 10e-3])

% Get depletion approximation solution
[JV_DA.tau1em7, r0_DA.tau1em7, Voc_DA.tau1em7, g0_DA.tau1em7] =  calcR0_pn_junction(par_pn_tau1em7);
