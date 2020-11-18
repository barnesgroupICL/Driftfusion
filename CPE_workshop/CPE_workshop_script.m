% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

%% Initialise the code. Only needs to be run once at the start of each session
initialise_df

%% Input parameters
params_filepath = './CPE_workshop_Input_files/blocking_contacts_3_layer.csv';   % Filepath to the parameters file
output_filename = 'blocking_contacts_3_layer';      % Filename for output file

light_intensities = [0, 1, 2, 3, 4, 5, 6];          % List the light intensities here
Vmax = 1.2;                                         % Maximum voltage for cyclic voltammogram
Vmin = -1.2;                                        % Minimum voltage for cyclic voltammogram

xlimits = [0, 0];   % Sets the limits for the CV plot x-axis. Set to [0 ,0] for autoscaling
ylimits = [0, 0];   % Sets the limits for the CV plot y-axis. Set to [0 ,0] for autoscaling

%% Call function to obtain equilibrium and cyclic voltammogram solutions
sol = sigma_light(params_filepath, light_intensities, Vmax, Vmin, output_filename);

%% Rename solution
eval([output_filename, '= sol;'])

%% Plot equilibrium band diagram
dfplot.ELx(sol.sol_eq.el)

%% Plot cyclic voltammograms
plot_CVs(sol.sol_CVs, light_intensities, xlimits, ylimits);

%% Plot conductivity (sigma) vs light intensity
plot_sigma_light(sol.sol_CVs, light_intensities, sol.sigma);