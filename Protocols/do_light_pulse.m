function [sol_pulse] = do_light_pulse(sol_ini, pulse_int, tmax, tpoints, duty, mobseti, log_timemesh)
% Uses square wave light generator for light source 2 and peforms a single
% pulse

%% Input arguments
% SOL_INI = initial conditions
% PULSE_INT
% MOBSETI = Ion mobility switch
% RS = Series resistance - recommended to use Rs = 1e6 for approx open
% circuit

disp('Starting light pulse')
par = sol_ini.par;

if log_timemesh
    par.tmesh_type = 2;
    par.tmax = tmax;
    par.t0 = tmax/1e4;
else
    % Setup time mesh
    par.tmesh_type = 1;
    par.tmax = tmax;
    par.t0 = 0;
end

par.mobseti = mobseti;

% Setup square wave function generator
par.g2_fun_type = 'square';
par.tpoints = tpoints;
par.g2_fun_arg(1) = 0;            % Lower intensity
par.g2_fun_arg(2) = pulse_int;    % Higher intensity
par.g2_fun_arg(3) = tmax+(1e-3*tmax);         % Capture length
par.g2_fun_arg(4) = duty;         % Duty cycle [%]

sol_pulse = df(sol_ini, par);

disp('light pulse complete')

end