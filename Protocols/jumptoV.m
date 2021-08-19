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
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%% Start code
if tdwell == 0
    error('tdwell cannot be set to zero- please choose a time step > 0 s')
end

disp('Starting jump to V')
%% OUTPUTS
% SOL_RELAX - the output relaxation solution
par = sol_ini.par;

%% Set up sweep
% Characteristic diffusion time
t_diff = (par.dcum0(end)^2)/(2*par.kB*par.T*min(min(par.mu_n), min(par.mu_p)));
par.tmax = 100*t_diff;
par.tmesh_type = 1;
par.t0 = 0;
par.V_fun_type = 'sweep';
par.V_fun_arg(1) = par.Vapp;
par.V_fun_arg(2) = Vjump;
par.V_fun_arg(3) = 100*t_diff;

par.mobseti = 0;

disp('Initial jump with zero ion mobility')
jump1 = df(sol_ini, par);
disp('Complete')

par = jump1.par;

if accelerate
    % Take ratio of electron and ion mobilities in the active layer
    rat_anion = par.mu_n(par.active_layer)/par.mu_a(par.active_layer);
    rat_cation = par.mu_n(par.active_layer)/par.mu_c(par.active_layer);
    
    % If the ratio is infinity (ion mobility set to zero) then set the ratio to
    % zero instead
    if isnan(rat_anion) || isinf(rat_anion)
        rat_anion = 0;
    end
    
    if isnan(rat_cation) || isinf(rat_cation)
        rat_cation = 0;
    end
    par.K_a = rat_anion;
    par.K_c = rat_cation;
else
    par.K_a = 1;
    par.K_c = 1;
end

par.V_fun_type = 'constant';
par.V_fun_arg(1) = par.Vapp;    % For future proof

par.g1_fun_type = 'constant';
par.g1_fun_arg(1) = Int;        % For future proof
par.int1 = Int;

par.mobseti = mobseti;
par.tmax = tdwell;
par.t0 = par.tmax/1e8;
par.tmesh_type = 2;
par.tpoints = 200;

disp('Dwell stage...')
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
disp('Complete')
sol_relax = sol;
sol_relax.par.K_c = 1;
sol_relax.par.K_a = 1;
% Read out currents from LH side
dfplot.Jt(sol_relax, sol_relax.x(end));

disp('Jump to V complete')
end