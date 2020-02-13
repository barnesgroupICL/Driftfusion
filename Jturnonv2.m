function [sol_turnon] = Jturnonv2(sol_sweep, tmax, tpoints, Vapp, Int)

par = sol_sweep.par;


par.mobseti = 1;
par.JV = 0;
par.Vapp = Vapp;
par.tmax = tmax;
par.tpoints = tpoints;

sol_turnon = df(sol_sweep, par);

end