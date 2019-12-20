% An example script for comparing the current-voltage curves for two different
% architects of perovskite solar cell
initialise_df

%% Create parameters objects
par_1a = pc('Input_files/3_layer_methods_test1a.csv');
par_1b = pc('Input_files/3_layer_methods_test1b.csv');
par_2a = pc('Input_files/3_layer_methods_test2a.csv');
par_2b = pc('Input_files/3_layer_methods_test2b.csv');

%% Find equilibrium solutions
soleq_1a = equilibrate(par_1a, 1);
soleq_1b = equilibrate(par_1b, 1);
soleq_2a = equilibrate(par_2a, 1);
soleq_2b = equilibrate(par_2b, 1);

%% Perform current-voltage scans
% JV = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_1a = doJV(soleq_1a.el, 1e-9, 261, 1, 0, 0, 1.3, 3);
JV_1b = doJV(soleq_1b.el, 1e-9, 261, 1, 0, 0, 1.3, 3);
JV_2a = doJV(soleq_2a.el, 1e-9, 261, 1, 0, 0, 1.3, 3);
JV_2b = doJV(soleq_2b.el, 1e-9, 261, 1, 0, 0, 1.3, 3);

%% Plot the current voltage curves
dfplot.JV(JV_1a,3)
hold on
dfplot.JV(JV_1b,3)
hold on
dfplot.JV(JV_2a,3)
hold on
dfplot.JV(JV_2b,3)
hold on
legend('PS 1a, dark, F', 'PS 1a, dark, R', 'PS 1a, 1 sun, F', 'PS 1a, 1 sun, R',...
    'PS 1b, dark, F', 'PS 1b, dark, R', 'PS 1b, 1 sun, F', 'PS 1b, 1 sun, R',...
    'PS 2a, dark, F', 'PS 2a, dark, R', 'PS 2a, 1 sun, F', 'PS 2a, 1 sun, R',...
    'PS 2b, dark, F', 'PS 2b, dark, R', 'PS 2b, 1 sun, F', 'PS 2b, 1 sun, R')
xlim([0, 1.3])
ylim([-25e-3, 10e-3])
hold off

%% Plot the electrostatic potential for PS 1a during forward illuminated scan
dfplot.Vx(JV_1a.ill.f, 1e9*[0, 0.2, 0.4, 0.6, 0.8]);
legend('HTL', 'interface 1', 'Active layer', 'Interface 2', 'ETL',...
    '0 V', '0.2 V', '0.4 V', '0.6 V', '0.8 V')

%% Plot the electron and hole profiles for PS 1a during forward illuminated scan
dfplot.npx(JV_1a.ill.f, 1e9*[0, 0.2, 0.4, 0.6, 0.8]);
ylim([1e9, 1e18])
legend('HTL','interface 1','Active layer','Interface 2','ETL',...
    'n, 0 V', 'p, 0 V', 'n, 0.2 V', 'p, 0.2 V', 'n, 0.4 V', 'p, 0.4 V',...
    'n, 0.6 V', 'p, 0.6 V', 'n, 0.8 V', 'p, 0.8 V')
