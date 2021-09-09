function sol_IMPS = doIMPS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
% Performs IMPS measurement at the
% specified FREQUENCY
% INT_BASE = constant (DC) bias component [mulitples of GX]
% INT_DELTA = variable (AC) bias component [mulitples of GX]
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
disp(['Starting IMPS, base intensity ', num2str(int_base), ', delta intensity ', num2str(int_delta),...
    ', at frequency ', num2str(frequency), ' Hz']);
if int_base < int_delta
    error('Constant bias intensity cannot be less than the AC intensity magnitude')
end
par = sol_ini.par;

% Setup time mesh
par.tmesh_type = 1;
par.tmax = tmax;
par.t0 = 0;

% Setup square wave function generator
par.g1_fun_type = 'sin';
par.tpoints = tpoints;
par.g1_fun_arg(1) = int_base;
par.g1_fun_arg(2) = int_delta;
par.g1_fun_arg(3) = frequency;
par.g1_fun_arg(4) = 0;          % phase

disp('Applying oscillating optical bias')
sol_IMPS = df(sol_ini, par);

disp('IMPS complete')
end


