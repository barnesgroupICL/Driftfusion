function sol = Vapp_function(sol_ini, Vapp_func, coeff, tmax, tpoints, logtime)

par = sol_ini.par;
par.Vapp_func = Vapp_func;
par.tmax = tmax;
par.tpoints = tpoints;
par.JV = 2;
par.Vapp_params = coeff;

if logtime
    par.tmesh_type = 2;
    par.t0 = par.tmax/1e6;
else
    par.tmesh_type = 1;
    par.t0 = 0;
end

sol = df(sol_ini, par);

end