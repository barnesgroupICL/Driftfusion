function sol_Jtran = doJtrans(sol_ini, tmax, pulselen, pulseint)

p = sol_ini.p;
p.tmesh_type = 2;
p.tmax = tmax;
p.t0 = p.tmax/1e6;
p.pulseon = 1;
p.pulseint = pulseint;
p.pulselen = pulselen;
p.pulsestart = 1e-8;

sol_Jtran = pindrift(sol_ini, p);
sol_Jtran.p.pulseon = 0;

end