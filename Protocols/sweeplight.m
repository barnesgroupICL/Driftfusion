function sol_sweep = sweeplight(sol_ini, tmax, tpoints, start_int, end_int)

par = sol_ini.par;
par.g1_fun_type = 'sweep';
par.tmesh_type = 2;
par.tmax = tmax;
par.t0 = tmax/1e8;
par.tpoints = tpoints;
par.gen_arg1(1) = start_int;
par.gen_arg1(2) = end_int;

sol_sweep = df(sol_ini, par);

end