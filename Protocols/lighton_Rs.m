function sol_ill = lighton_Rs(sol_ini, Int, stab_time, mobseti, Rs, pnts)
% A function to switch the light on and record the transient state with a
% series resistance Rs
%% Input arguments
% SOL_INI = initial conditions
% INT = Light intensity (Suns)
% STAB_TIME = Stabilisation time - length of the transient. Setting to -1
% to cycle to stable solution from initial time of 0.1
% MOBSETI = Ion mobility switch
% RS = Series resistance - recommended to use Rs = 1e6 for approx open
% circuit
% PNTS = Number of time points in the solution

disp(['Starting LIGHTON, Intensity = ', num2str(Int), ' Rs = ', num2str(Rs)])
par = sol_ini.par;
par_origin = par;

par.tmesh_type = 2;
par.tmax = 1e-3;
par.t0 = par.tmax/1e6;
par.tpoints = 10;
par.JV = 0;
par.mobseti = 0;
par.Rs = Rs;
par.Rs_initial = 1;

disp('Switching on series resistance- initial fast linear sweep')
sol_Rs = df(sol_ini, par);

disp('Switching on illumination')
par.Rs_initial = 0;
par.Int = Int;

par.mobseti = mobseti;
if stab_time >= 0
    par.tmax = stab_time;
else
    par.tmax = 0.1;
end
par.t0 = par.tmax/1e6;
par.tpoints = pnts;

sol = df(sol_Rs, par);

if stab_time == -1
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
end

sol_ill = sol;
sol_ill.par.mobseti = par_origin.mobseti;

%dfplot.ELx(sol_ill);
%dfplot.Voct(sol_ill);
disp('LIGHTON complete')

end