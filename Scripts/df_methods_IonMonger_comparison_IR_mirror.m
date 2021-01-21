initialise_df

%% Read in base parameters
par_df_im_mirror = pc('Input_files/IonMonger_default_IR_mirror.csv');

%% Get equilibrium solutions
soleq_df_im_mirror= equilibrate(par_df_im_mirror);

%% Perform JV scan at 1.2 Vs-1
% For comparison with IonMonger default parameters
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
CV_df_im_50mVs_mirror = doCV(soleq_df_im_mirror.ion, 1.0398, 0, 1.2, 0, 50e-3, 1, 241);
CV_df_im_100mVs_mirror = doCV(soleq_df_im_mirror.ion, 1.0398, 0, 1.2, 0, 100e-3, 1, 241);
CV_df_im_200mVs_mirror = doCV(soleq_df_im_mirror.ion, 1.0398, 0, 1.2, 0, 200e-3, 1, 241);

%% CV plots
dfplot.JtotVapp(CV_df_im_50mVs_mirror, 0)
hold on
dfplot.JtotVapp(CV_df_im_100mVs_mirror, 0)
hold on
dfplot.JtotVapp(CV_df_im_200mVs_mirror, 0)
ylim([-25e-3,10e-3])
xlim([0, 1.1])
legend('50 mVs-1', '100 mVs-1', '200 mVs-1') 
hold off