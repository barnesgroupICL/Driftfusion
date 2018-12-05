function sol_Jtran = doJtransient(sol_ini, pulseint, pulselen)
% Quick function to perform a J transient using the built-in pulse
% functionality

par = sol_ini.par;
par.pulseon = 1;
par.pulselen = pulselen;       % Transient pulse length
par.pulseint = pulseint;            % Transient pulse intensity- for BL and TM models, 100 mW/cm2 assumed
par.pulsestart = 0;%1e-6;
par.tmesh_type = 1;
par.tmax = 10e-6;
par.t0 = 0;%tmax/1000;
par.tpoints = 400;
par.calcJ = 0;

sol_Jtran = df(sol_ini, par);
sol_Jtran.par.pulseon = 0;

end