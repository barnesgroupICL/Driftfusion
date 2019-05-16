function sol_relax = jumptoV(sol_ini, Vjump, tdwell, Int)
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
tjump = 1e-6;   % Jump to V in 1 us

% Calculate linear rate coefficient
coeff = (V1 - V0)/tjump;
Vapp_func = @(coeff, t) coeff(1)*t;

sol_ini.par.mobseti = 0;
% Vapp_function(sol_ini, Vapp_func, tmax, tpoints, logtime)
jump1 = Vapp_function(sol_ini, Vapp_func, coeff, tjump, 20, 0);

par.Vapp = V1;
par.mobseti = 1;
par.tmax = tdwell;
par.t0 = 0;%par.tmax/1e6;
par.tmesh_type = 1;
par.tpoints = 200;
par.Int = Int;

sol_relax = df(jump1, par);

% Read out currents from LH side
dfplot.Jt(sol_relax, 1);

end