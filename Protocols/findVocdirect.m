function sol_Voc = findVocdirect(sol_ini, light_intensity)

par = sol_ini.par;

par.Rs = 1e6;
par.tmesh_type = 1;
par.t0 = 0;
par.tmax = 1e-6;
par.tpoints = 40;

sol_Rs = df(sol_ini, par);

par.int1 = light_intensity;
par.tmax = 1e-3;

sol_Voc = df(sol_Rs, par);

end