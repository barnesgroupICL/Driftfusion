function sol_relax = jumptoV(sol_ini, Vjump, tdwell, Int, K_cation, stabilise)
% A function to simulate a jump-to-Voltage measurement as used for SDP
%% Input arugments
% SOL_INI - input solution- could be cell at equilibrium or stablised at an
% applied voltage
% VJUMP - the voltage to be jumped to
% TDWELL - the relaxation time after jumping

%% OUTPUTS
% SOL_RELAX - the output relaxation solution
par = sol_ini.par;
V0 = par.Vapp;
V1 = Vjump;
tjump = 1e-10;   % Jump to V in 1 us

% Calculate linear rate coefficient
coeff = (V1 - V0)/tjump;
Vapp_func = @(coeff, t) coeff(1)*t;

sol_ini.par.mobseti = 0;
% Vapp_function(sol_ini, Vapp_func, tmax, tpoints, logtime)
jump1 = Vapp_function(sol_ini, Vapp_func, coeff, tjump, 20, 1);

par.K_cation = K_cation;
par.Vapp = V1;
par.mobseti = 1;
par.tmax = tdwell;
par.t0 = par.tmax/1e8;
par.tmesh_type = 2;
par.tpoints = 100;
par.Int = Int;

sol = df(jump1, par);

if stabilise
    all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    % loop to check ions have reached stable config- if not accelerate ions by
    while any(all_stable) == 0
        disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10^j)]);
        
        par.tmax = par.tmax*10;
        par.t0 = par.tmax/1e6;
        
        sol = df(sol, par);
        
        all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    end
end

sol_relax = sol;

sol_relax.par.K_cation = 1;
% Read out currents from LH side
dfplot.Jt(sol_relax, 1);

end