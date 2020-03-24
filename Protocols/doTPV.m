function [sol_TPV, sol_ill] = doTPV(sol_ini, bias_int, stab_time, mobseti, Rs, pulse_int, tmax, tpoints, duty)
%% Input arguments
% SOL_INI = initial conditions
% INT1 = Bias light intensity (Suns)
% STAB_TIME = Stabilisation time - length of the transient. Setting to -1
% to cycle to stable solution from initial time of 0.1
% MOBSETI = Ion mobility switch
% RS = Series resistance - recommended to use Rs = 1e6 for approx open
% circuit
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Start code
% Get open circuit solution
sol_ill = lightonRs(sol_ini, bias_int, stab_time, mobseti, Rs, tpoints);
% A function to switch the light on and record the transient state with a
% series resistance Rs
par = sol_ill.par;

% Setup time mesh
par.tmesh_type = 1;
par.t0 = 0;
par.tmax = tmax;

% Setup square wave function generator
par.g2_fun_type = 'square';
par.tpoints = tpoints;
par.g2_fun_arg(1) = 0;
par.g2_fun_arg(2) = pulse_int;
par.g2_fun_arg(3) = tmax;
par.g2_fun_arg(4) = duty;           % Duty cycle

sol_TPV = df(sol_ill, par);

end