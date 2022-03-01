%% Conductivity profiles

%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('Input_files/alox_2_gates.csv');
%% Equilibrium solutions
soleq_alox = equilibrate(par_alox);

%% Plot equilibrium energy level diagram
dfplot.ELnpx(soleq_alox.ion)

%% Current-voltage scan
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
% JVsol = doJV(soleq_sio2.ion, 100e-3, 201, 1, 0, 0, 1, 1);

% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
k_scan = 0.001;
Vmin = -15;
Vmax = 15;
sol_CV = doCV(soleq_alox.ion, 0, 0, Vmax, Vmin, k_scan, 1, 201);
%% Plot Vapp vs time
% dfplot.Vappt(sol_CV)

%% Plot JV scan
dfplot.JtotVapp(sol_CV, 0);
%set(gca,'YScale','log')

%% Plot anion and cation densities
%dfplot.acx(sol_CV, 1/k_scan*[0, 0.5, 1.0, 2.5, 3.0]);

%% Plot electron and hole profiles
dfplot.npx(sol_CV, 1/k_scan*[0, 0.5, 1.0, 2.5, 3.0]);

%% Plot space charge density
dfplot.rhox(sol_CV, 1/k_scan*[0, 0.5, 1.0, 2.5, 3.0]);

%% Calculate conductivity
[sigma_n, sigma_p] = dfana.calc_conductivity(sol_CV);

%% Debye length Calculation
par = par_alox;

e = par.e;
V_T = par.kB*par.T;                     % Thermal votlage
epp_pvsk = e*par.epp0*par.epp(3);       % Perovskite absolute dielectric constant
N0 = par.Ncat(3);                   
N0_courtier = 1.6e19;                   % cm-3

L_D = sqrt((epp_pvsk*V_T)/(e*N0));      % Deby width [cm]
L_D_courtier = sqrt((epp_pvsk*V_T)/(e*N0_courtier));

N_Debye = 3;                            % Number of Debye lengths to average electron density over
%%
x_perov_left = dcum0(3);

x = sol_CV.x;
t = sol_CV.t;
Vappt = dfana.calcVapp(sol_CV);
% Get point at which perovskite starts 
sigma_n_bar = mean(sigma_n(:, x > x_perov_left & x < x_perov_left + N_Debye*L_D), 2);
sigma_p_bar = mean(sigma_p(:, x > x_perov_left & x < x_perov_left + N_Debye*L_D), 2);

sigma_n_bar_entire = mean(sigma_n(:, x > x_perov_left & x < x_perov_left + 4.00E-05), 2);
sigma_p_bar_entire = mean(sigma_p(:, x > x_perov_left & x < x_perov_left + 4.00E-05), 2);

sigma_n_bar_bulk = mean(sigma_n(:, x > x_perov_left + N_Debye*L_D & x < x_perov_left + 4.00E-05), 2);
sigma_p_bar_bulk = mean(sigma_p(:, x > x_perov_left + N_Debye*L_D & x < x_perov_left + 4.00E-05), 2);



%% Find peak conductivity for applied bias
pp_Vmax = find(Vappt == max(Vappt));      %% pp = point position
pp_Vmin = find(Vappt == min(Vappt));      %% pp = point position

sigma_n_bar_Vpeak = sigma_n_bar(pp_Vmax);
sigma_p_bar_Vpeak = sigma_p_bar(pp_Vmax);

% %% Plot average conductivity
% figure
% plot(Vappt, sigma_n_bar, Vappt, sigma_p_bar)
% axis([-1 1 0 inf])
% xlabel('Voltage [V]')
% ylabel('Average conductivity [Siemens]')
% legend('Electron', 'Hole')

%%
%% Plot average conductivity
figure
semilogy(Vappt, sigma_n_bar, Vappt, sigma_p_bar)
xlabel('Voltage [V]')
ylabel('Average channel conductivity [Semilog]')
legend('Electron', 'Hole')

%% Plot average conductivity
figure
plot(Vappt, sigma_n_bar, Vappt, sigma_p_bar)
xlabel('Voltage [V]')
ylabel('Average channel conductivity [Linear]')
legend('Electron', 'Hole')

%%
%% Plot average conductivity
figure
semilogy(Vappt, sigma_n_bar_bulk, Vappt, sigma_p_bar_bulk)
xlabel('Voltage [V]')
ylabel('Average bulk conductivity [Semilog]')
legend('Electron', 'Hole')

%% Plot average conductivity
figure
plot(Vappt, sigma_n_bar_bulk, Vappt, sigma_p_bar_bulk)
xlabel('Voltage [V]')
ylabel('Average bulk conductivity [Linear]')
legend('Electron', 'Hole')

%% Plot average conductivity
figure
semilogy(Vappt, sigma_n_bar_entire, Vappt, sigma_p_bar_entire)
xlabel('Voltage [V]')
ylabel('Average entire conductivity [Semilog]')
legend('Electron', 'Hole')

%% Plot average conductivity
figure
plot(Vappt, sigma_n_bar_entire, Vappt, sigma_p_bar_entire)
xlabel('Voltage [V]')
ylabel('Average entire conductivity [Linear]')
legend('Electron', 'Hole')
%% Plot Peak conductivity
% PC - how do you intend to plot this? The peak voltage only occurs
% once per voltage cycle so you cannot plot as a function of voltage as you
% have tried below. 


% plot(Vappt, sigma_n_peak_positive_voltage, Vappt, sigma_p_peak_positive_voltage)
% xlabel('Voltage [V]')
% ylabel('Peak conductivity [Siemens]')
% legend('Electron', 'Hole')
%% Make movie for anions and cations
%makemovie(sol_CV, @dfplot.acx, 0, [0, 1.5e18], 'acx', true, true);