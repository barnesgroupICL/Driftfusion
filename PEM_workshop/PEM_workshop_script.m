% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

%% Initialise the code. Only needs to be run once at the start of each session
initialise_df

%% Input parameters
params_filepath = './PEM_workshop_Input_files/intrinsic_ohmic.csv';     % Filepath to the parameters file
output_filename = 'intrinsic_ohmic';                                    % Filename for output file

light_intensities = 0;          % List the light intensities here
Vmax = 1.2;                     % Maximum voltage for cyclic voltammogram
Vmin = -1.2;                    % Minimum voltage for cyclic voltammogram

xlimits = [0, 0];   % Sets the limits for the CV plot x-axis. Set to [0 ,0] for autoscaling
ylimits = [0, 0];   % Sets the limits for the CV plot y-axis. Set to [0 ,0] for autoscaling

%% Load in parameters
par = pc(params_filepath);

%% Obtain and plot equilibrium solution
sol_eq = equilibrate(par, 1);
% Plot equilibrium band diagram
dfplot.ELx(sol_eq.el);
% Plot carrier densities at equilbrium
dfplot.npx(sol_eq.el);
% Export equilibrium band daigram
export_solution(['Equilibrium_solution_', output_filename], sol_eq.el, 0)
export_ELx(['Equilibrium_energy_levels_', output_filename], sol_eq.el, 0)

%% Call function to obtain equilibrium and cyclic voltammogram solutions
sol_CVs = get_CVs(sol_eq.el, light_intensities, Vmax, Vmin, output_filename);

%% Output current and voltage for first light intensity (usually 0)
Vapp = dfana.calcVapp(sol_CVs.sol(1));  % Applied voltage
Jstruct = dfana.calcJ(sol_CVs.sol(1));
J = Jstruct.tot(:,1);                   % Total current at left-hand boundary

%% Rename solution
eval(['sol_eq_', output_filename, '= sol_eq;'])
eval(['sol_CVs_',output_filename, '= sol_CVs;'])

%% Export the solution for the first light intensity at specified voltage
Vapp = 0.6;
export_solution(['Solution at Vapp =', num2str(Vapp), output_filename], sol_CVs.sol(1), Vapp)

%% Plot cyclic voltammograms
plot_CVs(sol_CVs, light_intensities, xlimits, ylimits);

%% Plot conductivity (sigma) vs light intensity
plot_sigma_light(sol_CVs, light_intensities);