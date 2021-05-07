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
par_tmax = par;
attempted_frozen_ions = false;
% check if the solver was "Unable to meet integration
%  tolerances without reducing the step size below the smallest
%  value allowed at time t." and try to re-run the simulation with slightly
%  different parameters
while size(sol.u,1) ~= par.tpoints
    % as first attempt to fix, increase the tolerances
    if sol.par.RelTol < 1e-3
        par.RelTol = min(sol.par.RelTol*3.5, 1e-3);
        par.AbsTol = min(sol.par.AbsTol*3.5, 1e-5);
        warning('Driftfusion:VappFunction',...
            ['VappFunction the solver did not succeed, loosening the relative tolerance from '...
            num2str(sol.par.RelTol) ' to ' num2str(par.RelTol) ' and absolute tolerance from '...
            num2str(sol.par.AbsTol) ' to ' num2str(par.AbsTol)...
            '. BEWARE that the returned solution will have a high value of RelTol and AbsTol!'])
    % start the full simulation from a starting point with no moving ions
    % this changes the starting point, so it makes sense only in
    % simulations where the starting point is not important, like when
    % applying periodic potentials
    elseif ~attempted_frozen_ions && strcmp(Vapp_func, 'sin')
        attempted_frozen_ions = true;
        par_ions = par;
        par_ions.mobseti = false;
        warning('Driftfusion:VappFunction',...
            'VappFunction the solver did not succeed, trying to start from a simulation with frozen ions');
        sol_temp = df(sol_ini, par_ions);
        if size(sol_temp.u,1) == par_ions.tpoints
            sol_ini = sol_temp;
        end
    % start the full simulation from a short simulation
    % this changes the starting point, so it makes sense only in
    % simulations where the starting point is not important, like when
    % applying periodic potentials
    elseif i < 3 && strcmp(Vapp_func, 'sin')
        i = i + 1;
        old_tmax = par_tmax.tmax;
        par_tmax.tmax = old_tmax / 4;
        par_tmax.tpoints = round(par_tmax.tpoints / 4);
        warning('Driftfusion:VappFunction',...
            ['VappFunction the solver did not succeed, trying to start from a short simulation reducing tmax from '...
            num2str(old_tmax) ' to ' num2str(par_tmax.tmax)]);
        sol_temp = df(sol_ini, par_tmax);
        if size(sol_temp.u,1) == par_tmax.tpoints
            sol_ini = sol_temp;
        end
    % try to force smaller time steps
    elseif sol.par.MaxStepFactor > 0.001
        par.MaxStepFactor = sol.par.MaxStepFactor/10;
        warning('Driftfusion:VappFunction',...
            ['VappFunction the solver did not succeed, decreasing the maximum time step changing MaxStepFactor from '...
            num2str(sol.par.MaxStepFactor) ' to ' num2str(par.MaxStepFactor)...
            '. BEWARE that the returned solution will have a small value of MaxStepFactor!'])
    % desperately increase the tolerance
    elseif sol.par.RelTol < 1e-2
        par.RelTol = min(sol.par.RelTol*3.5, 1e-2);
        par.AbsTol = min(sol.par.AbsTol*3.5, 1e-4);
        warning('Driftfusion:VappFunction',...
            ['VappFunction the solver did not succeed, loosening the relative tolerance from '...
            num2str(sol.par.RelTol) ' to ' num2str(par.RelTol) ' and absolute tolerance from '...
            num2str(sol.par.AbsTol) ' to ' num2str(par.AbsTol)...
            '. BEWARE that the returned solution will have a high value of RelTol and AbsTol!'])
    else
        warning('Driftfusion:VappFunction',...
            ['VappFunction cannot make the simulation with parameters ' num2str(Vapp_coeff) ' work in any way'])
        break
    end
    % beware that the returned solution can have a value of RelTol, AbsTol
    % and MaxStepFactor different from the initially requested ones!
    sol = df(sol_ini, par);
end

end
