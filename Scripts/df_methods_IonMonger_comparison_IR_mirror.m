initialise_df

%% Read in base parameters
par_df_im_mirror= pc('Input_files/IonMonger_default_IR_mirror.csv');

%% Get equilibrium solutions
soleq_df_im_mirror= equilibrate(par_df_im);

%% Perform JV scan at 1.2 Vs-1
% For comparison with IonMonger default parameters
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
CV_df_im_50mVs = doCV(soleq_df_im.ion, 1.0398, 0, 1.2, 0, 50e-3, 1, 241);
CV_df_im_100mVs = doCV(soleq_df_im.ion, 1.0398, 0, 1.2, 0, 100e-3, 1, 241);
CV_df_im_200mVs = doCV(soleq_df_im.ion, 1.0398, 0, 1.2, 0, 200e-3, 1, 241);

%% CV plots
dfplot.JtotVapp(CV_df_im_50mVs, 0)
hold on
dfplot.JtotVapp(CV_df_im_100mVs, 0)
hold on
dfplot.JtotVapp(CV_df_im_200mVs, 0)
ylim([-10e-3,25e-3])
xlim([0, 1.1])
legend('50 mVs-1', '100 mVs-1', '200 mVs-1') 
hold off