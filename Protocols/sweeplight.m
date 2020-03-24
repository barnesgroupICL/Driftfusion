function sol_sweep = sweepLight(sol_ini, tmax, tpoints, end_int)
% Linear light intensity sweep from the intensity of the input solution SOL_INI to
% the intensity defined by END_INT
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

disp('Starting SWEEPLIGHT...')
par = sol_ini.par;
par.g1_fun_type = 'sweep';
par.tmesh_type = 2;
par.tmax = tmax;
par.t0 = tmax/1e8;
par.tpoints = tpoints;
par.g1_fun_arg(1) = par.int1;       % Use intensity from input sol
par.g1_fun_arg(2) = end_int;
par.g1_fun_arg(3) = tmax;

sol_sweep = df(sol_ini, par);
disp('Complete')
end