function sol_relax = jumptoV(sol_ini, Vjump, tdwell, mobseti, Int, stabilise, accelerate)
% A function to simulate a jump-to-Voltage measurement as used for SDP
%% Input arugments
% SOL_INI - input solution- could be cell at equilibrium or stablised at an
% applied voltage
% VJUMP - the voltage to be jumped to
% TDWELL - the relaxation time after jumping
% MOBSETI - determines whether the ion mobility is switched on after the
% jump
% INT - Intensity during stabilisation phase
% STABILISE - Determines whether the time step is increased such that the
% device is at steady-state at the end of the solution
% ACCELERATE - Accelerate the ions to be at the same mobility as the
% electrons (for easy access to steady state solutions)

if tdwell == 0
    error('tdwell cannot be set to zero- please choose a time step > 0 s')
end

disp('Starting jump to V')
%% OUTPUTS
% SOL_RELAX - the output relaxation solution
par = sol_ini.par;
V0 = par.Vapp;
V1 = Vjump;
par.tmesh_type = 1;
tjump = 1e-10;   % Jump to V in 1 us

% Calculate linear rate coefficient
coeff = [V0, ((V1 - V0)/tjump)];
Vapp_func = @(coeff, t) coeff(1) + coeff(2)*t;

sol_ini.par.mobseti = 0;
disp('Initial jump with zero ion mobility')
% Vapp_function(sol_ini, Vapp_func, tmax, tpoints, logtime)
jump1 = Vapp_function(sol_ini, Vapp_func, coeff, tjump, 20, 1);

if accelerate
    % Take ratio of electron and ion mobilities in the active layer
    rat_anion = par.mue(par.active_layer)/par.muion(par.active_layer);
    rat_cation = par.mue(par.active_layer)/par.mucat(par.active_layer);
    
    % If the ratio is infinity (ion mobility set to zero) then set the ratio to
    % zero instead
    if isnan(rat_anion) || isinf(rat_anion)
        rat_anion = 0;
    end
    
    if isnan(rat_cation) || isinf(rat_cation)
        rat_cation = 0;
    end
    par.K_anion = rat_anion;
    par.K_cation = rat_cation;
else
    par.K_anion = 1;
    par.K_cation = 1;
end

par.Vapp = V1;
par.mobseti = mobseti;
par.tmax = tdwell;
par.t0 = par.tmax/1e8;
par.tmesh_type = 2;
par.tpoints = 200;
par.int1 = Int;

sol = df(jump1, par);
j = 1;
if stabilise
    all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    % loop to check ions have reached stable config- if not accelerate ions by
    while any(all_stable) == 0       
        par.tmax = par.tmax*10;
        
        disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10^j)]);
        
        par.t0 = par.tmax/1e6;
        
        sol = df(jump1, par);
        
        all_stable = verifyStabilization(sol.u, sol.t, 0.7);
        j = j+1;
    end
end

sol_relax = sol;
sol_relax.par.K_cation = 1;
sol_relax.par.K_anion = 1;
% Read out currents from LH side
dfplot.Jt(sol_relax, sol_relax.x(end));

disp('Jump to V complete')
end