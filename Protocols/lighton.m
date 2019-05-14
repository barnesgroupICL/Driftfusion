function sol_ill = lighton(sol_ini, Int, stab_time, mobseti, Rs, pnts)
% A function to switch the light on and record the transient state with a
% series resistance Rs

disp(['Starting LIGHTON, Intensity = ', num2str(Int), ' Rs = ', num2str(Rs)])
par = sol_ini.par;
par_origin = par;

par.tmesh_type = 2;
par.tmax = 1e-3;
par.t0 = par.tmax/1e6;
par.tpoints = 40;
par.JV = 0;
par.mobseti = 0;
par.Rs = Rs;

disp('Switching on series resistance')
sol_Rs = df(sol_ini, par);

disp('Switching on illumination')
par.Int = Int;

par.mobseti = mobseti;
par.tmax = stab_time;
par.t0 = par.tmax/1e6;
par.tpoints = pnts;

sol = df(sol_Rs, par);

all_stable = verifyStabilization(sol.u, sol.t, 0.7);

% loop to check ions have reached stable config- if not accelerate ions by
% order of mag
j = 1;

while any(all_stable) == 0
    disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10^j)]);
    
    par.tmax = par.tmax*10;
    par.t0 = par.tmax/1e6;

    sol = df(sol, par);
    
    all_stable = verifyStabilization(sol.u, sol.t, 0.7);

end

sol_ill = sol;
sol_ill.par.mobseti = par_origin.mobseti;

dfplot.ELx(sol_ill);

disp('LIGHTON complete')

end