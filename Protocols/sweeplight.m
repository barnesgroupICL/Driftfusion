function sol_sweep = sweepLight(sol_ini, tmax, tpoints, end_int)
% Linear light intensity sweep from the intensity of the input solution SOL_INI to
% the intensity defined by END_INT
% P. Calado, Imperial College London, 2019

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