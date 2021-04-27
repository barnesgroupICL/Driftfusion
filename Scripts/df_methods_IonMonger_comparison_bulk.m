initialise_df

%% Read in base parameters
par_df_im_bulk = pc('Input_files/IonMonger_default_bulk.csv');

%% Get equilibrium solutions
soleq_df_im_bulk = equilibrate(par_df_im_bulk);

%% Perform JV scan at 1.2 Vs-1
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_0p1Vs_im = doCV(soleq_df_im_bulk.ion, 1.0398, 0, -1.3, 0, 0.1, 1, 261);
sol_CV_1Vs_im = doCV(soleq_df_im_bulk.ion, 1.0398, 0, -1.3, 0, 1, 1, 261);
sol_CV_10Vs_im = doCV(soleq_df_im_bulk.ion, 1.0398, 0, -1.3, 0, 10, 1, 261);

%% Plot JV - plot inverted due to layer ordering
dfplot.JtotVapp(sol_CV_0p1Vs_im, 0);
hold on
dfplot.JtotVapp(sol_CV_1Vs_im, 0);
hold on
dfplot.JtotVapp(sol_CV_10Vs_im, 0);
ylim([-10e-3,25e-3])
legend('0.1 Vs-1', '1 Vs-1', '10 Vs-1') 

%% Profile plots
% Electrostatic potential during 1 Vs-1 forward scan
dfplot.Vx(sol_CV_1Vs_im, [0, 0.2, 0.4, 0.6, 0.8])
legend
legend('HTL','interface 1','Active layer','Interface 2','ETL','0 V', '0.2 V', '0.4 V', '0.6 V', '0.8 V')

% Electronic carrier densities during during 1 Vs-1 forward scan
dfplot.npx(sol_CV_1Vs_im, [0, 0.2, 0.4, 0.6, 0.8])
ylim([1e10, 1e19])
legend('HTL','interface 1','Active layer','Interface 2','ETL',...
    'n, 0 V', 'p, 0 V', 'n, 0.2 V', 'p, 0.2 V', 'n, 0.4 V', 'p, 0.4 V',...
    'n, 0.6 V', 'p, 0.6 V', 'n, 0.8 V', 'p, 0.8 V')

% Ionic carrier densities during during 1 Vs-1 forward scan
dfplot.acx(sol_CV_1Vs_im, [0, 0.2, 0.4, 0.6, 0.8])
xlim([497, 502])
legend('HTL','interface 1','Active layer','Interface 2','ETL',...
    'a, 0 V', 'a, 0 V', 'a, 0.2 V', 'a, 0.2 V', 'a, 0.4 V', 'a, 0.4 V',...
    'c, 0.6 V', 'c, 0.6 V', 'c, 0.8 V', 'c, 0.8 V')

% Equilibrium
dfplot.ELx(soleq_df_im_bulk.ion)

