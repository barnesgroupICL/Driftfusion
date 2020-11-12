function CPE_master_fun(params_filepath, output_filename, light_intensities, Vmax, Vmin)

% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

%% Obtain cyclic voltammogram (CV) solutions and conductivity vs voltage (sigma)
sigma_sol = sigma_light(params_filepath, light_intensities, Vmax, Vmin);

%% Export CV solutions to tab delimited text file
headers = cell(1, length(light_intensities) + 1);
headers{1,1} = "Voltage";
for i=1:length(light_intensities)
    current_intensity ="Current, "+ num2str(light_intensities(i))+ " suns";
    headers{1,i+1} = current_intensity;
end

% Extract voltage and current arrays from the solution
V = dfana.calcVapp(sigma_sol.sol_CVs(1));
Js = zeros(length(V), length(light_intensities)); 
JV = zeros(length(V), length(light_intensities)+1);    % Eventually V and Js will  be transposed
JV(:,1) = V';   % Write the voltage array to first column

for i = 1:length(light_intensities)
   J = dfana.calcJ(sigma_sol.sol_CVs(i));     % J is structure containing all the current density elements
   JV(:, i+1) = J.tot(:, 1)';
end 

% Create the table and write to file
output_cell = [headers; num2cell(JV)];
T = cell2table(output_cell);
writetable(T, ['./Output_files/', output_filename, '.txt'], 'Delimiter', 'tab', 'WriteVariableNames', 0);

end