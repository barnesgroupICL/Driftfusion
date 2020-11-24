function [Vapp, R_A, sigma] = get_sigma(sol_CV)
% Created for the Centre for Processable Electronics (CPE) MRes course
% on semiconductor physics
% Authors: P.R.F. Barnes, P. Calado 2020
% Department of Physics, Imperial College London

% A function to calculate the conductivity of a CV
%% Input arguments
% SOL_CV is a solution from a cyclic voltammogram simulation
% SIGMA is a vector containing the conductivity calculated from the CV at each voltage
%% Outputs
% Vapp is a an array with the applied voltages
% R is an array with the inferred resistance at each voltage
% sigma is an array with the inferred conductivity at each voltage 

Vapp = dfana.calcVapp(sol_CV);
J = dfana.calcJ(sol_CV);            % [Acm-2] Calculates currents at left-hand boundary
R_A = gradient(Vapp, J.tot(:,1));   % [Ohms cm2] Area normalised resistance from dVdJ
d = sol_CV.par.dcum(end);           % [cm] Total device thickness
sigma = d./R_A;                     % [S cm-1] Conductivity 
end