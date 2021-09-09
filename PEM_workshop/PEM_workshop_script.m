% Semiconductor device physics computation simulation workshop
% Created for the Centre for Processable Electronics (CPE) MRes PEM course
% Department of Physics, Imperial College London
% Authors: P. R. F. Barnes, P. Calado
% Nov 2020

%% Input parameters
params_filepath = './PEM_workshop_Input_files/intrinsic_ohmic.csv';     % Filepath to the parameters file
output_filename = 'intrinsic_ohmic';                                    % Filename for output file
light_intensity = 1;            % Suns equivalent
Vmax = 1.2;                     % Maximum voltage for cyclic voltammogram
Vmin = -1.2;                    % Minimum voltage for cyclic voltammogram
scan_rate = 1e-3;               % Current-voltage scan rate [Vs-1] 

%% Load in parameters
par = pc(params_filepath);

%% Obtain and plot equilibrium solution
sol_eq = equilibrate(par, 1);

%% Call function to obtain equilibrium and cyclic voltammogram solutions
sol_CV = doCV(sol_eq.el, light_intensity, 0, Vmax, Vmin, scan_rate, 1, 401);

%% Output current and voltage for first light intensity (usually 0)
Vapp = dfana.calcVapp(sol_CV);  % Applied voltage
Jstruct = dfana.calcJ(sol_CV);
J = Jstruct.tot(:,1);                   % Total current at left-hand boundary

%% Plot cyclic voltammograms
xlimits = [0, 0];   % Sets the limits for the CV plot x-axis. Set to [0 ,0] for autoscaling
ylimits = [0, 0];   % Sets the limits for the CV plot y-axis. Set to [0 ,0] for autoscaling
plot_CVs(sol_CV, light_intensity, xlimits, ylimits);

%% Plot energy level diagram at applied bias Vapp for first light intensity
Vplot = 0;
% Get corresponding time, TPLOT for VPLOT
if Vplot >= 0
    tplot = Vplot/scan_rate;
elseif Vplot < 0
    tplot = ((2*Vmax)+abs(Vplot))/scan_rate; 
end
% PLot the energy level diagram at time TPLOT
dfplot.ELx(sol_CV, tplot);
% PLot the carrier densities at time TPLOT
dfplot.npx(sol_CV, tplot);

%% Rename solution
eval(['sol_eq_', output_filename, '= sol_eq;'])
eval(['sol_CV_',output_filename, '= sol_CV;'])

%% Exporting the solutions
% Export the solution for the first light intensity at voltage VPLOT
export_solution(['Solution at Vapp =', num2str(Vplot), ' V ', output_filename], sol_CV, tplot)
export_ELx(['Energy levels at Vapp =', num2str(Vplot), ' V ', output_filename], sol_CV, tplot)
