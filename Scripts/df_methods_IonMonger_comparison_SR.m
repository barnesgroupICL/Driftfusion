initialise_df

%% Read in base parameters
par_df_im = pc('Input_files/IonMonger_default_IR.csv');

%% Get equilibrium solutions
soleq_df_im = equilibrate(par_df_im);

%% Perform JV scan at 1.2 Vs-1
% For comparison with IonMonger default parameters
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
CV_df_50mVs_df = doCV(soleq_df_im.ion, 1.0398, 0, -1.2, 0, 50e-3, 1, 241);
CV_df_100mVs_df = doCV(soleq_df_im.ion, 1.0398, 0, -1.2, 0, 100e-3, 1, 241);
CV_df_200mVs_df = doCV(soleq_df_im.ion, 1.0398, 0, -1.2, 0, 200e-3, 1, 241);

%       doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
% JV_df_50mVs_im = doJV(soleq_df_im.ion, 50e-3, 121, 1.0398, 1, 0, -1.2, 2);
% JV_df_100mVs_im = doJV(soleq_df_im.ion, 100e-3, 121, 1.0398, 1, 0, -1.2, 2);
% JV_df_200mVs_im = doJV(soleq_df_im.ion, 200e-3, 121, 1.0398, 1, 0, -1.2, 2);

%% Plot JV - plot inverted due to layer ordering
% dfplot.JV(JV_df_50mVs_im, 2);
% hold on
% dfplot.JV(JV_df_100mVs_im, 2);
% hold on
% dfplot.JV(JV_df_200mVs_im, 2);
% legend('50 mVs-1, F', '50 mVs-1, R', '100 mVs-1, F', '100 mVs-1, R', '200 mVs-1, F', '200 mVs-1, R') 

dfplot.JtotVapp(CV_df_50mVs_df, 0)
hold on
dfplot.JtotVapp(CV_df_100mVs_df, 0)
hold on
dfplot.JtotVapp(CV_df_200mVs_df, 0)
ylim([-10e-3,25e-3])
xlim([0, 1.1])
legend('50 mVs-1', '100 mVs-1', '200 mVs-1') 
hold off

% %% Electrostatic potential during 1 Vs-1 forward scan
% dfplot.Vx(JV_df_1V_im.ill.f, [0, 0.2, 0.4, 0.6, 0.8])
% legend
% legend('HTL','interface 1','Active layer','Interface 2','ETL','0 V', '0.2 V', '0.4 V', '0.6 V', '0.8 V')
% 
% %% Electronic carrier densities during during 1 Vs-1 forward scan
% dfplot.npx(JV_df_1V_im.ill.f, [0, 0.2, 0.4, 0.6, 0.8])
% ylim([1e10, 1e19])
% legend('HTL','interface 1','Active layer','Interface 2','ETL',...
%     'n, 0 V', 'p, 0 V', 'n, 0.2 V', 'p, 0.2 V', 'n, 0.4 V', 'p, 0.4 V',...
%     'n, 0.6 V', 'p, 0.6 V', 'n, 0.8 V', 'p, 0.8 V')
% 
% %% Ionic carrier densities during during 1 Vs-1 forward scan
% dfplot.acx(JV_df_1V_im.ill.f, [0, 0.2, 0.4, 0.6, 0.8])
% xlim([497, 502])
% legend('HTL','interface 1','Active layer','Interface 2','ETL',...
%     'a, 0 V', 'a, 0 V', 'a, 0.2 V', 'a, 0.2 V', 'a, 0.4 V', 'a, 0.4 V',...
%     'c, 0.6 V', 'c, 0.6 V', 'c, 0.8 V', 'c, 0.8 V')
% 
% %% Equilibrium
% dfplot.ELx(soleq_df_im.ion)

