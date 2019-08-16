function sol_IMPS = do_IMPS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
% Switches to open circuit at INT_BASE, runs to stabilised OC then performs IMPS measurement at the
% specified FREQUENCY
% INT_BASE = constant (DC) bias component [mulitples of GX]
% INT_DELTA = variable (AC) bias component [mulitples of GX]
disp(['Starting IMPS, base intensity ', num2str(int_base), ', delta intensity ', num2str(int_delta),...
    ', at frequency ', num2str(frequency), ' Hz']);
par = sol_ini.par;

% Setup time mesh
par.tmesh_type = 1;
par.tmax = tmax;
par.t0 = 0;

% Setup square wave function generator
par.g1_fun_type = 'sin';
par.tpoints = tpoints;
par.g1_fun_arg(1) = int_base;
par.g1_fun_arg(2) = int_delta;
par.g1_fun_arg(3) = frequency;
par.g1_fun_arg(4) = 0;          % phase

disp('Applying oscillating optical bias')
sol_IMPS = df(sol_ini, par);

%% Plot the J as a function of t
dfplot.Jt(sol_IMPS,0)

disp('IMPS complete')
end


