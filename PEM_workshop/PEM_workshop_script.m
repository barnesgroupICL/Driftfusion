% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

%% Initialise the code. Only needs to be run once at the start of each session
initialise_df

%% Input parameters
params_filepath = './PEM_workshop_Input_files/intrinsic_ohmic.csv';     % Filepath to the parameters file
output_filename = 'intrinsic_ohmic';                                    % Filename for output file
light_intensities = [0];        % List the light intensities here in square brackets e.g. [0, 1, 2]
Vmax = 1.2;                     % Maximum voltage for cyclic voltammogram
Vmin = -1.2;                    % Minimum voltage for cyclic voltammogram
scan_rate = 1e-3;               % Current-voltage scan rate [Vs-1] 
xlimits = [0, 0];   % Sets the limits for the CV plot x-axis. Set to [0 ,0] for autoscaling
ylimits = [0, 0];   % Sets the limits for the CV plot y-axis. Set to [0 ,0] for autoscaling

%% Load in parameters
par = pc(params_filepath);

%% Obtain and plot equilibrium solution
sol_eq = equilibrate(par, 1);

%% Call function to obtain equilibrium and cyclic voltammogram solutions
sol_CV = get_CVs(sol_eq.el, light_intensities, Vmax, Vmin, scan_rate, output_filename);

%% Output current and voltage for first light intensity (usually 0)
Vapp = dfana.calcVapp(sol_CV.sol(1));  % Applied voltage
Jstruct = dfana.calcJ(sol_CV.sol(1));
J = Jstruct.tot(:,1);                   % Total current at left-hand boundary

%% Rename solution
eval(['sol_eq_', output_filename, '= sol_eq;'])
eval(['sol_CV_',output_filename, '= sol_CV;'])

%% Plot cyclic voltammograms
plot_CVs(sol_CV, light_intensities, xlimits, ylimits);

%% Plot energy level diagram at applied bias Vapp for first light intensity
Vplot = 0;
% Get corresponding time, TPLOT for VPLOT
if Vplot >= 0
    tplot = Vplot/scan_rate;
elseif Vplot < 0
    tplot = ((2*Vmax)+abs(Vplot))/scan_rate; 
end
% PLot the energy level diagram at time TPLOT
dfplot.ELx(sol_CV.sol(1), tplot);
% PLot the carrier densities at time TPLOT
dfplot.npx(sol_CV.sol(1), tplot);

%% Exporting the solutions
% Export equilibrium electron, hole and electrostatic potential densities
export_solution(['Equilibrium_solution_', output_filename], sol_eq.el, 0)

% Export the solution for the first light intensity at specified voltage
Vexport = 0.6;
export_solution(['Solution at Vapp =', num2str(Vexport), ' V ', output_filename], sol_CV.sol(1), Vexport)
