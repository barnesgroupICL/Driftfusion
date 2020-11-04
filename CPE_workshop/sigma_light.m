function [sol_eq, sol_CVs, light_intensities, sigma] = sigma_light(params_filepath, light_intensities, Vmax, Vmin)
% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

%% Input parameters
% PARAMS_FILEPATH

%% Load in parameters
par = pc(params_filepath);

%% Obtain equilibrium solution
sol_eq = equilibrate(par);

%% Cyclic-voltammogram initial parameters
V0 = 0;         % Start scan voltage [V]
Vmax = 2;       % Max scan voltage [V]
Vmin = -2;      % min scan voltage [V]
scan_rate = 1e-3;  % [Vs-1]
cycles = 1;
Vpoints = 401;  % Number of voltage points in output solution

%% preallocate memory
sigma = zeros(length(light_intensities), Vpoints);
R = zeros(length(light_intensities), Vpoints);

for i = 1:length(light_intensities)
    % Run cyclic voltammogram
    sol_CVs(i) = doCV(sol_eq.el, light_intensities(i), V0, Vmax, Vmin, scan_rate, cycles, Vpoints);
    
    % Calculate resistance and conductivity
    [Vapp(i,:), R(i,:), sigma(i,:)] = get_sigma(sol_CVs(i));
end

end

