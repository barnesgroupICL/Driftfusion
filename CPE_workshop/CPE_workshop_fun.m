function [sol_eq, sol_CVs, light_intensities, sigma] = get_sigmas(params_filepath, light_intensities)
% Test script for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

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
sigma = zeros(lenght(light_intensities), Vpoints);
R = zeros(lenght(light_intensities), Vpoints);

for i = 1:length(light_intensities)
    % Run cyclic voltammogram
    sol_CVs(i) = doCV(soleq_intrinsic_ohmic.el, light_intensities(i), V0, Vmax, Vmin, scan_rate, cycles, Vpoints);
    
    % Calculate resistance and conductivity
    [Vapp(i,:), R(i,:), sigma(i,:)] = get_sigma(sol_CVs(i));
end

end

