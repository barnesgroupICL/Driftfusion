function [Vapp, R, sigma] = get_sigma(sol_CV)
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
J = dfana.calcJ(sol_CV);        % Calculates currents at left-hand boundary
R = gradient(Vapp, J.tot(:,1));             % Area normalised resistance from J-V characteristic
sigma = sol_CV.par.dcum(end)./R;    % sigma = d/R
end