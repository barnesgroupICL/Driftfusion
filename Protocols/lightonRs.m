function sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
% A function to switch the light on and record the transient state with a
% series resistance Rs
%% Input arguments
% SOL_INI = initial conditions
% int1 = Bias light intensity (Suns)
% stable_time = Stabilisation time - length of the transient. Setting to -stable_time
% to cycle to stable solution from initial time of stable_time
% MOBSETI = Ion mobility switch
% RS = Series resistance - recommended to use Rs = 1e6 for approx open
% circuit
% PNTS = Number of time points in the solution
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
disp(['Starting LIGHTON_RS, Intensity = ', num2str(int1), ' Rs = ', num2str(Rs)])
par = sol_ini.par;
par_origin = par;

disp('Switching on illumination')
par.mobseti = 0;
par.tmax = 1e-3;
par.t0 = par.tmax/1e6;
par.Rs_initial = 0;
par.int1 = int1;
par.g1_fun_type = 'sweep';
% COEFF = [Amplitude_initial, Amplitude_final, tmax]
par.g1_fun_arg = [0, int1, par.tmax];
par.JV = 0;

sol_ill1 = df(sol_ini, par);

par.g1_fun_type = 'constant';

% Characteristic diffusion time
t_diff = (par.dcum0(end)^2)/(2*par.kB*par.T*min(min(par.mue), min(par.muh)));

C = par.e*par.epp(par.active_layer)*par.epp0/(par.d_active);      % Parallel plate capacitance
tau_RC = C*Rs;

par.tmesh_type = 2;
par.tmax = tau_RC;
par.t0 = par.tmax/1e6;
par.tpoints = 20;

if Rs > 0
    par.Rs = Rs;
    par.Rs_initial = 1;
    
    disp('Switching on series resistance- initial linear sweep')
    sol = df(sol_ill1, par);
    disp('Complete')
    
    % Longer linear time step to reach stabilisation
    par.Rs_initial = 0;
    par.tmesh_type = 2;
    par.tmax = tau_RC*1e-3;
    par.t0 = par.tmax/1e6;
    
    disp('Stabilisation step with Rs')
    sol_Rs = df(sol, par);
    disp('Complete')
else
    sol_Rs = sol_ill1;
end

par.mobseti = mobseti;
% If STABLE_TIME is entered as zero set to default value
if stable_time == 0
    par.tmax = 10*t_diff;
else
	par.tmax = abs(stable_time);
end
par.t0 = par.tmax/1e6;
par.tpoints = pnts;

disp(['Performing initial transient with tmax = ', num2str(par.tmax)]);
sol = df(sol_Rs, par);

% Switch off precondition voltage
par.Vapp = 0;

% Stable time is enetered as negative then run to stable solution
if stable_time < 0
    sol = stabilize(sol);
end

sol_ill = sol;
sol_ill.par.mobseti = par_origin.mobseti;

disp('LIGHTON_RS complete')

end
