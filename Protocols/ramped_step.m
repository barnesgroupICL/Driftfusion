function [sol_ramp, sol_dwell] = ramped_step(sol_in, deltaV, t_ramp, t_dwell)

par = sol_in.par;

% Set up time mesh
par.tmesh_type = 1;
par.t0 = 0;
par.tmax = t_ramp;
par.tpoints = 100;

%% Define the voltage function
% Sweep coefficients
% COEFF = [Amplitude_initial, Amplitude_final, tmax]
par.V_fun_type = 'sweep';           
par.V_fun_arg(1) = par.Vapp;           
par.V_fun_arg(2) = par.Vapp + deltaV;      
par.V_fun_arg(3) = t_ramp;                   

disp('Applying ramped potential')
sol_ramp = df(sol_in, par);

par = sol_ramp.par;

par.tmesh_type = 2;
par.tmax = t_dwell;
par.t0 = t_dwell/1e8;
par.tpoints = 300;

par.V_fun_type = 'constant';           
par.V_fun_arg(1) = par.Vapp;   

disp('Applying constant potential')
sol_dwell = df(sol_ramp, par);

% Concatenate solutions
sol_concat = sol_dwell;
sol_concat.u = [sol_ramp.u; sol_dwell.u];
sol_concat.t = [sol_ramp.t, sol_ramp.t(end) + sol_dwell.t];

end