% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

%% Intrinsic device with Ohmic contacts
params_filepath = './Input_files/intrinsic_ohmic.csv';
light_intensities = [0, 1, 2, 3, 4, 5, 6];   % List the light intensities here
Vmax = 1;       % Maximum voltage for cyclic voltammogram
Vmin = -1;      % Minimum voltage for cyclic voltammogram

%% Get equilibrium (eq) and cyclic voltammogram (CV) solutions
[sol_eq, sol_CVs, light_intensities, sigma] = sigma_light(params_filepath, light_intensities, Vmax, Vmin);

%% plot CVs
plot_CVs(sol_CVs, light_intensities)

%% plot conductivity (sigma) vs light intensity
plot_sigma_light(sol_CVs, light_intensities, sigma);