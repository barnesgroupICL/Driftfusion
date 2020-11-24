% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

%% Initialise the code. Only needs to be run once at the start of each session
initialise_df

%% Input parameters
params_filepath = './CPE_workshop_Input_files/blocking_contacts_3_layer.csv';   % Filepath to the parameters file
output_filename = 'blocking_contacts_3_layer';      % Filename for output file

light_intensities = [0, 1, 2, 3];          % List the light intensities here
Vmax = 1.3;                                         % Maximum voltage for cyclic voltammogram
Vmin = -1.3;                                        % Minimum voltage for cyclic voltammogram

xlimits = [0, 0];   % Sets the limits for the CV plot x-axis. Set to [0 ,0] for autoscaling
ylimits = [0, 0];   % Sets the limits for the CV plot y-axis. Set to [0 ,0] for autoscaling

%% Load in parameters
par = pc(params_filepath);

%% Obtain and plot equilibrium solution
sol_eq = equilibrate(par, 1);

% Plot equilibrium band diagram
dfplot.ELx(sol_eq.el)

% Export equilibrium band daigram
export_ELx(['Equilibrium_energy_levels_', output_filename], sol_eq.el, sol_eq.el.t(end))

%% Call function to obtain equilibrium and cyclic voltammogram solutions
sol_CVs = get_CVs(sol_eq.el, light_intensities, Vmax, Vmin, output_filename);

%% Rename solution
eval(['sol_eq_', output_filename, '= sol_eq;'])
eval(['sol_CVs_',output_filename, '= sol_CVs;'])

%% Plot cyclic voltammograms
plot_CVs(sol_CVs, light_intensities, xlimits, ylimits);

%% Plot conductivity (sigma) vs light intensity
plot_sigma_light(sol_CVs, light_intensities);