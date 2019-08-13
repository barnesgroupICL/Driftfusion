function sol_Jtran = doJtransient(sol_ini, tmax, tperiod, tpoints, pulseint)
% Quick function to perform a J transient using the built-in pulse
% functionality

par = sol_ini.par;
par.tmesh_type = 1;
par.tmax = tmax;
par.t0 = tmax/1e6;

par.g2_fun_type = 'square';
par.tpoints = tpoints;
par.gen_arg2(1) = 0;
par.gen_arg2(2) = pulseint;
par.gen_arg2(3) = tperiod;
par.gen_arg2(4) = 50;           % Duty cycle

sol_Jtran = df(sol_ini, par);
sol_Jtran.par.pulseon = 0;

dfplot.Jt(sol_Jtran, 1);

end