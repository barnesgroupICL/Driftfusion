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
Vmax = 10;       % Max scan voltage [V]
Vmin = -10;      % min scan voltage [V]
scan_rate = 1;  % [Vs-1]
cycles = 1;
Vpoints = 201;  % Number of voltage points in output solution

%% Run cyclic voltammogram
sol_0sun = doCV(soleq_intrinsic_ohmic.el, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, Vpoints);

%% Increase light intensity and run cyclic voltammogram
light_intensity = 1; % [Suns]
sol_1sun = doCV(soleq_intrinsic_ohmic.el, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, Vpoints);

%% plot solution
dfplot.JtotVapp(sol_0sun, 0)      % Second argument is the position at which the current is calculated - in this example it is the left-hand boundary
hold on
dfplot.JtotVapp(sol_1sun, 0)

    