function sol_IMVS = do_IMVS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
% Switches to open circuit at INT_BASE, runs to stabilised OC then performs IMVS measurement at the
% specified FREQUENCY
% INT_BASE = constant (DC) bias component [mulitples of GX]
% INT_DELTA = variable (AC) bias component [mulitples of GX]
disp(['Starting IMVS, base intensity ', num2str(int_base), ', delta intensity ', num2str(int_delta),...
    ', at frequency ', num2str(frequency)]);
par = sol_ini.par;

%sol_ill = lighton_Rs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
sol_ill = lighton_Rs(sol_ini, int_base, -1, 1, 1e6, 100);
par = sol_ill.par;

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
sol_IMVS = df(sol_ill, par);

%% Plot the Voc as a function of t
dfplot.Voct(sol_IMVS);

disp('IMVS complete')

end


