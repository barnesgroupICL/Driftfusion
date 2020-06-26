function sol_TPV = doTPVonly(sol_ill, pulse_int, tmax, tpoints, duty, tmesh_type)
%% Input arguments
% SOL_INI = initial conditions
% INT1 = Bias light intensity (Suns)
% STAB_TIME = Stabilisation time - length of the transient. Setting to -1
% to cycle to stable solution from initial time of 0.1
% MOBSETI = Ion mobility switch
% RS = Series resistance - recommended to use Rs = 1e6 for approx open
% circuit
% PULSE_INT = Pulse light intensity
% TMAX = Total length of the capture including pulse and decay
% DUTY = Duty cycle: Percentage of TMAX for which the pulse light is on
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

par = sol_ill.par;
% Setup time mesh
par.tmesh_type = tmesh_type;
switch tmesh_type
    case 1
        assert(1/tpoints < duty/100, 'Not even two time points while the laser pulse is on! increase tpoints or increase duty or use tmesh_type 2.')
        par.t0 = 0;
    case 2
        par.t0 = tmax*duty/100/10;
end
par.tmax = tmax;

% Setup square wave function generator
par.g2_fun_type = 'square';
par.tpoints = tpoints;
par.g2_fun_arg(1) = 0;
par.g2_fun_arg(2) = pulse_int;
% the square function defined in fun_gen adds an "on" point at the end of
% arg(3), which is not wanted in TPV, so I make arg(3) slightly larger than
% tmax
par.g2_fun_arg(3) = tmax*(1+eps);
par.g2_fun_arg(4) = duty;           % Duty cycle

sol_TPV = df(sol_ill, par);

end