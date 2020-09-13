function sol = VappFunction(sol_ini, Vapp_func, Vapp_coeff, tmax, tpoints, logtime)
% Applies a voltage function using the initial conditions defined by
% SOL_INI
% Input arguments:
% SOL_INI   = Solution containing the intial conditions
% VAPP_FUN  = 'constant', 'sweep', 'square', 'sin'
% VAPP_COEFF = coefficients array- see FUN_GEN for coefficients definition
% for each function type
% TMAX      = Final time position
% TPOINTS   = Number of solution time points (must be >=2 points
% LOGTIME   = 0 = linear time mesh, 1 = log time mesh.
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
par = sol_ini.par;
par.V_fun_type = Vapp_func;
par.tmax = tmax;
par.tpoints = tpoints;
par.V_fun_arg = Vapp_coeff;

if logtime
    par.tmesh_type = 2;
    par.t0 = par.tmax/1e6;
else
    par.tmesh_type = 1;
    par.t0 = 0;
end

sol = df(sol_ini, par);

i = 0;
while size(sol.u,1) ~= par.tpoints
    warning('Driftfusion:VappFunction',...
        'VappFunction the solver did not succeed, increasing time points')
    par.tpoints = par.tpoints * 2;
    sol = df(sol_ini, par);
    i = i+1;
    % in case 10 repetitions were not enough, stop
    if i > 4
        warning('Driftfusion:VappFunction',...
            'VappFunction increased the time points but the simulation broke anyway')
        break
    end
end

end