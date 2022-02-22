%% Conductivity profiles

%% Initialize driftfusion
initialise_df

%% Add parameter file to path 
% Filepath Mac
par_alox = pc('Input_files/alox.csv');
%% Equilibrium solutions
soleq_alox = equilibrate(par_alox);

%% Plot equilibrium energy level diagram
dfplot.ELnpx(soleq_alox.ion)

%% Current-voltage scan
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
% JVsol = doJV(soleq_sio2.ion, 100e-3, 201, 1, 0, 0, 1, 1);

% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
k_scan = 0.001;
Vmax=1;
sol_CV = doCV(soleq_alox.ion, 0, 0, Vmax, 0, k_scan, 1, 201);
%% Plot Vapp vs time
dfplot.Vappt(sol_CV)


%% Plot JV scan
dfplot.JtotVapp(sol_CV, 0);
%set(gca,'YScale','log')

%% Plot anion and cation densities
dfplot.acx(sol_CV, 1/k_scan*[0, 0.5, 1.0, 2.5, 3.0]);

%% Plot electron and hole profiles
dfplot.npx(sol_CV, 1/k_scan*[0, 0.5, 1.0, 2.5, 3.0]);

%% Plot space charge density
dfplot.rhox(sol_CV, 1/k_scan*[0, 0.5, 1.0, 2.5, 3.0]);

%% Calculate conductivity
[sigma_n, sigma_p] = dfana.calc_conductivity(sol_CV);
[sigma_n_bar_peak_positive_voltage, sigma_p_bar_peak_positive_voltage] = dfana.calc_peak_conductivity(sol_CV);

%% Debye length Calculation
L_D = 30e-7;
Thermal_voltage=25.7;
e0=8.854E-12; %in cm
Permittivity_perovskite=24.1.*e0;
N0=par_alox.Ncat;
N0_ownvalue=1e25;
e=par_alox.e;
Debye_Length = sqrt((Permittivity_perovskite.*25.7)./(e.*N0));
Debye_Length_own = sqrt((Permittivity_perovskite.*25.7)./(e.*N0_ownvalue));
%%
x_perov_left = 202e-7;
x = sol_CV.x;
t = sol_CV.t;
Vappt = dfana.calcVapp(sol_CV);
% Get point at which perovskite starts 
sigma_n_bar = mean(sigma_n(:, x > x_perov_left & x < x_perov_left + N_Debye*L_D), 2);
sigma_p_bar = mean(sigma_p(:, x > x_perov_left & x < x_perov_left + N_Debye*L_D), 2);

%% Find peak conductivity for applied bias
sigma_n_peak_positive_voltage = mean( sigma_n_bar_peak_positive_voltage(:, x > x_perov_left & x < x_perov_left + N_Debye*L_D), 2);
sigma_p_peak_positive_voltage = mean( sigma_p_bar_peak_positive_voltage(:, x > x_perov_left & x < x_perov_left + N_Debye*L_D), 2);

%% Plot average conductivity
figure
plot(Vappt, sigma_n_bar, Vappt, sigma_p_bar)
xlabel('Voltage [V]')
ylabel('Average conductivity [Siemens]')
legend('Electron', 'Hole')

 %% Plot Peak conductivity

plot(Vappt, sigma_n_peak_positive_voltage, Vappt, sigma_p_peak_positive_voltage)
xlabel('Voltage [V]')
ylabel('Peak conductivity [Siemens]')
legend('Electron', 'Hole')
%% Make movie for anions and cations
%makemovie(sol_CV, @dfplot.acx, 0, [0, 1.5e18], 'acx', true, true);