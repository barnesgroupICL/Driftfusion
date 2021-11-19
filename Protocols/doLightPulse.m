function [sol_pulse] = doLightPulse(sol_ini, pulse_int, tmax, tpoints, duty, mobseti, log_timemesh)
% Uses square wave light generator for light source 2 and peforms a single
% pulse
%
%% Input arguments
% SOL_INI = initial conditions
% PULSE_INT
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
disp('Starting light pulse')
par = sol_ini.par;

if log_timemesh
    par.tmesh_type = 2;
    par.tmax = tmax;
    par.t0 = tmax/1e4;
else
    % Setup time mesh
    par.tmesh_type = 1;
    par.tmax = tmax;
    par.t0 = 0;
end

par.mobseti = mobseti;

%% Setup square wave function generator
par.g2_fun_type = 'square';
par.tpoints = tpoints;
par.g2_fun_arg(1) = 0;                      % Lower intensity
par.g2_fun_arg(2) = pulse_int;              % Higher intensity
par.g2_fun_arg(3) = tmax+(1e-3*tmax);       % Capture length
par.g2_fun_arg(4) = duty;                   % Duty cycle [%]

sol_pulse = df(sol_ini, par);

disp('light pulse complete')

end