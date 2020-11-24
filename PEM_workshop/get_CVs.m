function sol_CVs = get_CVs(sol_CVs_in, light_intensities, Vmax, Vmin, output_filename)
% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London
%% Input arguments
% VMAX = max scan voltage [V]
% VMIN = min scan voltage [V]

%% Cyclic-voltammogram initial parameters
V0 = 0;         % Start scan voltage [V]
scan_rate = 1e-3;  % [Vs-1]
cycles = 1;
Vpoints = 401;  % Number of voltage points in output sol_CVsution

%% preallocate memory
sol_CVs.sigma = zeros(length(light_intensities), Vpoints);
sol_CVs.R = zeros(length(light_intensities), Vpoints);

for i = 1:length(light_intensities)
    % Run cyclic voltammogram
    sol_CVs.sol(i) = doCV(sol_CVs_in, light_intensities(i), V0, Vmax, Vmin, scan_rate, cycles, Vpoints);
    
    % Calculate resistance and conductivity
    [sol_CVs.Vapp(i,:), sol_CVs.R(i,:), sol_CVs.sigma(i,:)] = get_sigma(sol_CVs.sol(i));
end

% Export the CVs to .txt file in the Output_files folder
export_CVs(output_filename, light_intensities, sol_CVs.sol)

end

