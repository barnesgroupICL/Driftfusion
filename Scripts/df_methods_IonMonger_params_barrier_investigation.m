initialise_df

%% Read in base parameters
par_df_im_60p = pc('Input_files/IonMonger_default_IR_60p.csv');

%% Get equilibrium solutions
soleq_df_im_60p = equilibrate(par_df_im_60p);

%% Perform JV scan at 1.2 Vs-1
% For comparison with IonMonger default parameters
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
% CV_df_50mVs_im_60p = doCV(soleq_df_im_60p.ion, 1.0398, 0, -1.2, 0, 50e-3, 1, 241);
% CV_df_100mVs_im_60p = doCV(soleq_df_im_60p.ion, 1.0398, 0, -1.2, 0, 100e-3, 1, 241);
CV_df_200mVs_im_60p = doCV(soleq_df_im_60p.ion, 1.0398, 0, -1.2, 0, 200e-3, 1, 241);

%% Read in base parameters
par_df_im_100p = pc('Input_files/IonMonger_default_IR_100p.csv');

%% Get equilibrium solutions
soleq_df_im_100p = equilibrate(par_df_im_100p);

%% Perform JV scan at 1.2 Vs-1
% For comparison with IonMonger default parameters
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
% CV_df_50mVs_im_100p = doCV(soleq_df_im_100p.ion, 1.0398, 0, -1.2, 0, 50e-3, 1, 241);
% CV_df_100mVs_im_100p = doCV(soleq_df_im_100p.ion, 1.0398, 0, -1.2, 0, 100e-3, 1, 241);
CV_df_200mVs_im_100p = doCV(soleq_df_im_100p.ion, 1.0398, 0, -1.2, 0, 200e-3, 1, 241);

%% CV plots
% dfplot.JtotVapp(CV_df_50mVs_im_60p, 0)
% hold on
% dfplot.JtotVapp(CV_df_100mVs_im_60p, 0)
% hold on
dfplot.JtotVapp(CV_df_200mVs_im_60p, 0)
hold on

% dfplot.JtotVapp(CV_df_50mVs_im_100p, 0)
% hold on
% dfplot.JtotVapp(CV_df_100mVs_im_100p, 0)
% hold on
dfplot.JtotVapp(CV_df_200mVs_im_100p, 0)

figure(91)
hold on
% plotJV_im(sol_light_50mVs);
% plotJV_im(sol_light_100mVs);
plotJV_im(sol_light_200mVs);

ylim([-10e-3,25e-3])
xlim([-1.2, 0])

%% legend for comparison
legend('DF- 50 mVs-1, 60 pnt interface', 'DF- 100 mVs-1, 60 pnt interface', 'DF- 200 mVs-1,  60 pnt interface',...
   'DF- 50 mVs-1, 100 pnt interface', 'DF- 100 mVs-1, 100 pnt interface', 'DF- 200 mVs-1,  100 pnt interface',...
'IM- 50 mVs-1', 'IM- 100 mVs-1', 'IM- 200 mVs-1')