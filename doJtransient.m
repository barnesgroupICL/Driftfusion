function sol_Jtran = doJtransient(sol_ini)

p = sol_ini.p;
p.pulseon = 1;
p.pulselen = 10e-6;       % Transient pulse length
p.pulseint = 0.002;            % Transient pulse intensity- for BL and TM models, 100 mW/cm2 assumed
p.pulsestart = 0;%1e-6;
p.tmesh_type = 1;
p.tmax = 10e-6;
p.t0 = 0;%tmax/1000;
p.tpoints = 400;
p.deltat = p.tmax/(1e4*p.tpoints);
p.calcJ = 0;

sol_Jtran = df(sol_ini, p);
sol_Jtran.p.pulseon = 0;

end