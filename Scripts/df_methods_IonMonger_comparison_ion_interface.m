initialise_df

%% Read in base parameters
par_df_im = pc('Input_files/IonMonger_default_noIR.csv');

%% Get equilibrium solutions
soleq_df_im_0p85 = equilibrate(par_df_im);

%% Perform JV scan at 1.2 Vs-1
%       doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_df_0p1V = doJV(soleq_df.ion, 0.1, 131, 1.0398, 1, 0, -1.3, 2);
JV_df_1V_im = doJV(soleq_df_im.ion, 1, 131, 1.0398, 1, 0, -1.3, 2);
JV_df_10V = doJV(soleq_df.ion, 10, 131, 1.0398, 1, 0, -1.3, 2);

%% Plot JV
dfplot.JV(JV_df_1V_im, 2);
hold on
ylim([-10e-3,25e-3])

% %% Electrostatic potential during forward scan
% dfplot.Vx(JV_df_1V_im.ill.f, [JV_df_1V_im.ill.f.t(1), JV_df_1V_im.ill.f.t(21),...
%     JV_df_1V_im.ill.f.t(41), JV_df_1V_im.ill.f.t(61), JV_df_1V_im.ill.f.t(81)]);
% 
% %% Plot carrier densities during forward scan
% dfplot.npx(JV_df_1V_im.ill.f, [JV_df_1V_im.ill.f.t(1), JV_df_1V_im.ill.f.t(21),...
%     JV_df_1V_im.ill.f.t(41), JV_df_1V_im.ill.f.t(61), JV_df_1V_im.ill.f.t(81)]);
% ylim([1, 1e19])
% 
%% Plot ionic carrier densities during forward scan
dfplot.acx(JV_df_1V_im.ill.f, [JV_df_1V_im.ill.f.t(1), JV_df_1V_im.ill.f.t(21),...
    JV_df_1V_im.ill.f.t(41), JV_df_1V_im.ill.f.t(61), JV_df_1V_im.ill.f.t(81)]);
% 
% 
% %% Electrostatic potential during forward scan
% dfplot.Vx(JV_df_0p1V.ill.f, [JV_df_0p1V.ill.f.t(1), JV_df_0p1V.ill.f.t(21),...
%     JV_df_0p1V.ill.f.t(41), JV_df_0p1V.ill.f.t(61), JV_df_0p1V.ill.f.t(81)]);
% 
% %% Plot carrier densities during forward scan
% dfplot.npx(JV_df_0p1V.ill.f, [JV_df_0p1V.ill.f.t(1), JV_df_0p1V.ill.f.t(21),...
%     JV_df_0p1V.ill.f.t(41), JV_df_0p1V.ill.f.t(61), JV_df_0p1V.ill.f.t(81)]);
% ylim([1, 1e19])
% 
% %% Plot ionic carrier densities during forward scan
% dfplot.acx(JV_df_0p1V.ill.f, [JV_df_0p1V.ill.f.t(1), JV_df_0p1V.ill.f.t(21),...
%     JV_df_0p1V.ill.f.t(41), JV_df_0p1V.ill.f.t(61), JV_df_0p1V.ill.f.t(81)]);
% 
%% Equilibrium
dfplot.Vx(soleq_df_im.ion)
dfplot.npx(soleq_df_im.ion)
dfplot.acx(soleq_df_im.ion)

