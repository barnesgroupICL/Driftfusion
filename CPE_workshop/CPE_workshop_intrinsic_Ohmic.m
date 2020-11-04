% Test script for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London
%% Add folders and load defaults
initialise_df

%% Load in parameters
intrinsic_ohmic = pc('./Input_files/intrinsic_ohmic.csv');

%% Obtain equilibrium solution
soleq_intrinsic_ohmic = equilibrate(intrinsic_ohmic);

% Plot energy level diagram for the equilibrium solution
dfplot.ELx(soleq_intrinsic_ohmic.el)

%% Cyclic-voltammogram initial parameters
light_intensity = 0; % [Suns]
V0 = 0;         % Start scan voltage [V]
Vmax = 2;       % Max scan voltage [V]
Vmin = -2;      % min scan voltage [V]
scan_rate = 1;  % [Vs-1]
cycles = 1;
Vpoints = 401;  % Number of voltage points in output solution

%% Run cyclic voltammogram
intrinsic_ohmic_CV_0sun = doCV(soleq_intrinsic_ohmic.el, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, Vpoints);

%% Increase light intensity and run cyclic voltammogram
light_intensity = 1; % [Suns]
intrinsic_ohmic_CV_1sun = doCV(soleq_intrinsic_ohmic.el, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, Vpoints);

%% Increase light intensity and run cyclic voltammogram
light_intensity = 2; % [Suns]
intrinsic_ohmic_CV_2sun = doCV(soleq_intrinsic_ohmic.el, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, Vpoints);

%% plot solution
dfplot.JtotVapp(intrinsic_ohmic_CV_0sun, 0)      % Second argument is the position at which the current is calculated - in this example it is the left-hand boundary
hold on
dfplot.JtotVapp(intrinsic_ohmic_CV_1sun, 0)
dfplot.JtotVapp(intrinsic_ohmic_CV_2sun, 0)
legend('0 sun', '1 sun', '2 sun')

%% Calculate conductivity
[R_0sun, sigma_0sun] = get_sigma(intrinsic_ohmic_CV_0sun);
[R_1sun, sigma_1sun] = get_sigma(intrinsic_ohmic_CV_1sun);
[R_2sun, sigma_2sun] = get_sigma(intrinsic_ohmic_CV_2sun);

%%
Vapp = dfana.calcVapp(intrinsic_ohmic_CV_0sun);
figure(100)
plot(Vapp, sigma_0sun, Vapp, sigma_1sun, Vapp, sigma_2sun)
xlabel('Applied voltage [V]')
ylabel('Inferred Conductivity [S]')
legend('0 sun', '1 sun', '2 sun')

%%
figure(101)
plot(Vapp, R_0sun, Vapp, R_1sun, Vapp, R_2sun)
xlabel('Applied voltage [V]')
ylabel('Inferred Resistance [Ohms]')
legend('0 sun', '1 sun', '2 sun')
%% Make movie of carrier densities
%npx_intrinsic_0sun = makemovie(intrinsic_ohmic_CV_0sun, @dfplot.npx, 0, [1e10, 1e11], 'npx_intrinsic_0sun', 1, 0);
%npx_intrinsic_10sun = makemovie(intrinsic_ohmic_CV_1sun, @dfplot.npx, 0, 0, 'npx_intrinsic_10sun', 1, 0);
%Fx_intrinsic_0sun = makemovie(intrinsic_ohmic_CV_0sun, @dfplot.Fx, 0, 0, 'Fx_intrinsic_0sun', 1, 0);
%ELx_intrinsic_0sun = makemovie(intrinsic_ohmic_CV_0sun, @dfplot.ELx, 0, [-13, 7], 'ELx_intrinsic_0sun', 1, 0);


