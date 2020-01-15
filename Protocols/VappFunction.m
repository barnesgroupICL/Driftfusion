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

end