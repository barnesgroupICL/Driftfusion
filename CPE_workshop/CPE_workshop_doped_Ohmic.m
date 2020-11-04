% Test script for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London
%% Add folders and load defaults
initialise_df

%% Load in parameters
doped_ohmic = pc('./Input_files/doped_ohmic.csv');

%% Obtain equilibrium solution
soleq_doped_ohmic = equilibrate(doped_ohmic);

% Plot energy level diagram for the equilibrium solution
dfplot.ELx(soleq_doped_ohmic.el)

%% Cyclic-voltammogram initial parameters
light_intensity = 0; % [Suns]
V0 = 0;         % Start scan voltage [V]
Vmax = 2;       % Max scan voltage [V]
Vmin = -2;      % min scan voltage [V]
scan_rate = 1;  % [Vs-1]
cycles = 1;
Vpoints = 401;  % Number of voltage points in output solution

%% Run cyclic voltammogram
doped_ohmic_CV_0sun = doCV(soleq_doped_ohmic.el, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, Vpoints);

%% Increase light intensity and run cyclic voltammogram
light_intensity = 10; % [Suns]
doped_ohmic_CV_10sun = doCV(soleq_doped_ohmic.el, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, Vpoints);

%% Increase light intensity and run cyclic voltammogram
light_intensity = 20; % [Suns]
doped_ohmic_CV_20sun = doCV(soleq_doped_ohmic.el, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, Vpoints);

%% plot solution
dfplot.JtotVapp(doped_ohmic_CV_0sun, 0)      % Second argument is the position at which the current is calculated - in this example it is the left-hand boundary
hold on
dfplot.JtotVapp(doped_ohmic_CV_10sun, 0)
dfplot.JtotVapp(doped_ohmic_CV_20sun, 0)
legend('0 sun', '10 sun', '20 sun')
xlim([-0.1, 0.1])
hold off

%% Make movie of carrier densities
%npx_doped_0sun = makemovie(doped_ohmic_CV_0sun, @dfplot.npx, 0, [1e10, 1e11], 'npx_doped_0sun', 1, 0);
%npx_doped_10sun = makemovie(doped_ohmic_CV_10sun, @dfplot.npx, 0, 0, 'npx_doped_10sun', 1, 0);
%Fx_doped_0sun = makemovie(doped_ohmic_CV_0sun, @dfplot.Fx, 0, 0, 'Fx_doped_0sun', 1, 0);
%ELx_doped_0sun = makemovie(doped_ohmic_CV_0sun, @dfplot.ELx, 0, [-13, 7], 'ELx_doped_0sun', 1, 0);