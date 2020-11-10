% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

%% Input parameters
params_filepath = './CPE_workshop_Input_files/doped_ohmic.csv';  % Filepath to the appropriate parameters file
output_file_name = 'doped_ohmic';           % Filename for output file

light_intensities = [0, 1, 2, 3];      % List the light intensities here
Vmax = 2;                                       % Maximum voltage for cyclic voltammogram
Vmin = -2;                                      % Minimum voltage for cyclic voltammogram

%% Load in parameters
par = pc(params_filepath);

%% Obtain and plot equilibrium solution
sol_eq = equilibrate(par);

%% Obtain cyclic voltammogram (CV) solutions and conductivity vs voltage (sigma)
[sol_CVs, light_intensities, sigma] = sigma_light(sol_eq, light_intensities, Vmax, Vmin);

%% Export solutions to tab delimited text file
headers = cell(1, length(light_intensities) + 1);
headers{1,1} = "Voltage";
for i=1:length(light_intensities)
    current_intensity ="Current, "+ num2str(light_intensities(i))+ " suns";
    headers{1,i+1} = current_intensity;
end
%headers = cellstr(headers);

% Extract voltage and current arrays from the solution
V = dfana.calcVapp(sol_CVs(1));
Js = zeros(length(V), length(light_intensities)); 
JV = zeros(length(V), length(light_intensities)+1);    % Eventually V and Js will  be transposed
JV(:,1) = V';   % Write the voltage array to first column

for i = 1:length(light_intensities)
   J = dfana.calcJ(sol_CVs(i));     % J is structure containing all the current density elements
   JV(:, i+1) = J.tot(:, 1)';
end 

% Create the table and write to file
output_cell = [headers; num2cell(JV)];
T = cell2table(output_cell);
writetable(T, ['./Output_files/', output_file_name, '.txt'], 'Delimiter', 'tab', 'WriteVariableNames', 0);

%% Plot equilibrium band diagram
dfplot.ELx(sol_eq.el)

%% Plot cyclic voltammograms
xlimits = [-1, 1];      % Sets the limits for the plot x-axis
ylimits = [-0.1, 0.1];  % Sets the limits for the plot y-axis
plot_CVs(sol_CVs, light_intensities, xlimits, ylimits);

%% Plot conductivity (sigma) vs light intensity
plot_sigma_light(sol_CVs, light_intensities, sigma);