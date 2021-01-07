initialise_df

%% Read in base parameters
par_df_im = pc('Input_files/IonMonger_default_IR.csv');

%% Get equilibrium solutions
soleq_df_im = equilibrate(par_df_im);

%% JV method from Courtier, 2019, EES
sol_relax = jumptoV(soleq_df_im.ion, -1.2, 1e-3, 1, 0, 1, 1);
JV_im_sr = doJV(sol_relax, 100e-3, 241, 1.0398, 1, -1.2, 0, 2);

%% plot
dfplot.JV(JV_im_sr, 2)
xlim([-1.1,0])
ylim([-5e-3, 25e-3])
%% Perform JV scan at 1.2 Vs-1
%       doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JV_df_0p1V_im = doJV(soleq_df_im.ion, 0.1, 131, 1.0398, 1, 0, -1.3, 2);
JV_df_1V_im = doJV(soleq_df_im.ion, 1, 131, 1.0398, 1, 0, -1.3, 2);
JV_df_10V_im = doJV(soleq_df_im.ion, 10, 131, 1.0398, 1, 0, -1.3, 2);

%% Plot JV - plot inverted due to layer ordering
dfplot.JV(JV_df_0p1V_im, 2);
hold on
dfplot.JV(JV_df_1V_im, 2);
hold on
dfplot.JV(JV_df_10V_im, 2);
ylim([-10e-3,25e-3])
legend('0.1 Vs-1, F', '0.1 Vs-1, R', '1 Vs-1, F', '1 Vs-1, R', '10 Vs-1, F', '10 Vs-1, R') 

%% Electrostatic potential during 1 Vs-1 forward scan
dfplot.Vx(JV_df_1V_im.ill.f, [0, 0.2, 0.4, 0.6, 0.8])
legend
legend('HTL','interface 1','Active layer','Interface 2','ETL','0 V', '0.2 V', '0.4 V', '0.6 V', '0.8 V')

%% Electronic carrier densities during during 1 Vs-1 forward scan
dfplot.npx(JV_df_1V_im.ill.f, [0, 0.2, 0.4, 0.6, 0.8])
ylim([1e10, 1e19])
legend('HTL','interface 1','Active layer','Interface 2','ETL',...
    'n, 0 V', 'p, 0 V', 'n, 0.2 V', 'p, 0.2 V', 'n, 0.4 V', 'p, 0.4 V',...
    'n, 0.6 V', 'p, 0.6 V', 'n, 0.8 V', 'p, 0.8 V')

%% Ionic carrier densities during during 1 Vs-1 forward scan
dfplot.acx(JV_df_1V_im.ill.f, [0, 0.2, 0.4, 0.6, 0.8])
xlim([497, 502])
legend('HTL','interface 1','Active layer','Interface 2','ETL',...
    'a, 0 V', 'a, 0 V', 'a, 0.2 V', 'a, 0.2 V', 'a, 0.4 V', 'a, 0.4 V',...
    'c, 0.6 V', 'c, 0.6 V', 'c, 0.8 V', 'c, 0.8 V')

%% Equilibrium
dfplot.ELx(soleq_df_im.ion)

