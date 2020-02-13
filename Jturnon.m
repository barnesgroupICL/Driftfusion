function [Jturnonsweep, sol_turnon] = Jturnon(sol_ini, tmax, tpoints, Vapp, Int)

par = sol_ini.par;

% sweep to applied voltage with stationary ions @ 1 V/s
par.JV = 1;
Jturnonsweep = doJV(sol_ini, 1, 80, Int, 0, par.Vapp, Vapp, 2);

sol = Jturnonsweep.ill.f;

par.mobseti = 1;
par.JV = 0;
par.Vapp = Vapp;
par.tmax = tmax;
par.tpoints = tpoints;

sol_turnon = df(sol, par);

end