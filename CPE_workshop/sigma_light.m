function sigma_sol = sigma_light(params_filepath, light_intensities, Vmax, Vmin, output_filename)
% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London
%% Input arguments
% VMAX = max scan voltage [V]
% VMIN = min scan voltage [V]

%% Load in parameters
par = pc(params_filepath);

%% Obtain and plot equilibrium solution
sigma_sol.sol_eq = equilibrate(par);

%% Cyclic-voltammogram initial parameters
V0 = 0;         % Start scan voltage [V]
scan_rate = 1e-3;  % [Vs-1]
cycles = 1;
Vpoints = 401;  % Number of voltage points in output solution

%% preallocate memory
sigma_sol.sigma = zeros(length(light_intensities), Vpoints);
sigma_sol.R = zeros(length(light_intensities), Vpoints);

for i = 1:length(light_intensities)
    % Run cyclic voltammogram
    sigma_sol.sol_CVs(i) =...
        doCV(sigma_sol.sol_eq.el, light_intensities(i), V0, Vmax, Vmin, scan_rate, cycles, Vpoints);
    
    % Calculate resistance and conductivity
    [sigma_sol.Vapp(i,:), sigma_sol.R(i,:), sigma_sol.sigma(i,:)] =...
        get_sigma(sigma_sol.sol_CVs(i));
end

% Export the CVs to .txt file in the Output_files folder
export_CVs(output_filename, light_intensities, sigma_sol.sol_CVs)

end

