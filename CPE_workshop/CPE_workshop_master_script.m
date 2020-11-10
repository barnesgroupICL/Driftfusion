% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

%% Intrinsic device with Ohmic contacts
params_filepath = './Input_files/intrinsic_ohmic.csv';
light_intensities = [0, 1, 2, 3, 4, 5, 6];      % List the light intensities here
Vmax = 1;                                       % Maximum voltage for cyclic voltammogram
Vmin = -1;                                      % Minimum voltage for cyclic voltammogram
output_file_name = 'intrinsic_ohmic.txt';           % Filename for output file

%% Get equilibrium (eq) and cyclic voltammogram (CV) solutions
[sol_eq, sol_CVs, light_intensities, sigma] = sigma_light(params_filepath, light_intensities, Vmax, Vmin);

%% Export solutions to text file
headers = "Voltage";
for i=1:length(light_intensities)
    current_intensity ="Current, "+ num2str(light_intensities(i))+ " suns";
    headers = [headers, current_intensity];
end

%% Extract voltage and current arrays from the solution
V = dfana.calcVapp(sol_CVs(1));
Js = zeros(length(light_intensities), length(V)); 
JV = zeros(length(light_intensities)+1, length(V));    % Eventually V and Js will  be transposed
JV(1,:) = V;   % Write the voltage array to first column

for i = 1:length(light_intensities)
   J = dfana.calcJ(sol_CVs(i));     % J is structure containing all the current density elements
   JV(i+1,:) = J.tot(:,1);
end

T = table(JV, 'varnames', headers);

%writetable(headers, JV)

%% plot CVs
xlimits = [-0.1, 0.1];
plot_CVs(sol_CVs, light_intensities, xlimits);

%% plot conductivity (sigma) vs light intensity
plot_sigma_light(sol_CVs, light_intensities, sigma);