initialise_df

%% Read in base parameters
par_df_im = pc('Input_files/IonMonger_default_IR.csv');

%% Get equilibrium solutions
soleq_df_im = equilibrate(par_df_im);

%% Perform JV scan at 1.2 Vs-1
% For comparison with IonMonger default parameters
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
CV_df_50mVs_im = doCV(soleq_df_im.ion, 1.0398, 0, -1.2, 0, 50e-3, 1, 241);
CV_df_100mVs_im = doCV(soleq_df_im.ion, 1.0398, 0, -1.2, 0, 100e-3, 1, 241);
CV_df_200mVs_im = doCV(soleq_df_im.ion, 1.0398, 0, -1.2, 0, 200e-3, 1, 241);

%% CV plots
dfplot.JtotVapp(CV_df_50mVs_im, 0)
hold on
dfplot.JtotVapp(CV_df_100mVs_im, 0)
hold on
dfplot.JtotVapp(CV_df_200mVs_im, 0)
hold on
ylim([-10e-3,25e-3])
xlim([-1.2, 0])

%% Read in base parameters
par_df_im_mirror = pc('Input_files/IonMonger_default_IR_mirror.csv');

%% Get equilibrium solutions
soleq_df_im_mirror = equilibrate(par_df_im_mirror);

%% Perform JV scan at 1.2 Vs-1
% For comparison with IonMonger default parameters
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
CV_df_50mVs_im_mirror = doCV(soleq_df_im_mirror.ion, 1.0398, 0, 1.2, 0, 50e-3, 1, 241);
CV_df_100mVs_im_mirror = doCV(soleq_df_im_mirror.ion, 1.0398, 0, 1.2, 0, 100e-3, 1, 241);
CV_df_200mVs_im_mirror = doCV(soleq_df_im_mirror.ion, 1.0398, 0, 1.2, 0, 200e-3, 1, 241);

%% CV plots
dfplot.JtotVapp_mirror(CV_df_50mVs_im_mirror, 0)
hold on
dfplot.JtotVapp_mirror(CV_df_100mVs_im_mirror, 0)
hold on
dfplot.JtotVapp_mirror(CV_df_200mVs_im_mirror, 0)
ylim([-10e-3,25e-3])
xlim([-1.2, 0])
legend('50 mVs-1', '100 mVs-1', '200 mVs-1',...
    '50 mVs-1, mirror', '100 mVs-1, mirror', '200 mVs-1, mirror')
hold on

%% Profile plots
% Electrostatic potential during 1 Vs-1 forward scan
dfplot.Vx(CV_df_100mVs_im, 10*[0, 0.3, 0.6, 0.9])
legend('HTL','interface 1','Active layer','Interface 2','ETL','0 V', '0.3 V', '0.6 V', '0.9 V')

% Electronic carrier densities during during 1 Vs-1 forward scan
dfplot.npx(CV_df_100mVs_im, 10*[0, 0.3, 0.6, 0.9])
ylim([1e10, 1e19])
legend('HTL','interface 1','Active layer','Interface 2','ETL',...
    'n, 0 V', 'p, 0 V', 'n, 0.3 V', 'p, 0.3 V', 'n, 0.6 V', 'p, 0.6 V',...
    'n, 0.9 V', 'p, 0.9 V')

% Ionic carrier densities during during 1 Vs-1 forward scan
dfplot.acx(CV_df_100mVs_im, 10*[0, 0.3, 0.6, 0.9])
xlim([497, 502])
legend('HTL','interface 1','Active layer','Interface 2','ETL',...
    'a, 0 V', 'c, 0 V', 'Nani, 0 V', 'Ncat, 0 V', 'a, 0.3 V', 'c, 0.3 V', 'Nani, 0.3 V', 'Ncat, 0.3 V',...
    'a, 0.6 V', 'c, 0.6 V', 'Nani, 0.6 V', 'Ncat, 0.6 V', 'a, 0.9 V', 'c, 0.9 V', 'Nani, 0.9 V', 'Ncat, 0.9 V')

