function sol_ill = lighton_Rs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
% A function to switch the light on and record the transient state with a
% series resistance Rs
%% Input arguments
% SOL_INI = initial conditions
% int1 = Bias light intensity (Suns)
% stable_time = Stabilisation time - length of the transient. Setting to -1
% to cycle to stable solution from initial time of 0.1
% MOBSETI = Ion mobility switch
% RS = Series resistance - recommended to use Rs = 1e6 for approx open
% circuit
% PNTS = Number of time points in the solution

disp(['Starting LIGHTON, Intensity = ', num2str(int1), ' Rs = ', num2str(Rs)])
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

sol_ill1 = df(sol_ini, par);

par.g1_fun_type = 'constant';
par.tmesh_type = 2;
par.tmax = 1e-3;
par.t0 = par.tmax/1e6;
par.tpoints = 10;
par.JV = 0;

par.Rs = Rs;
par.Rs_initial = 1;

disp('Switching on series resistance- initial fast linear sweep')

sol_Rs = df(sol_ill1, par);

par.Rs_initial = 0;
par.mobseti = mobseti;
% If STABLE_TIME is  is entered as zero set to default value
if stable_time >= 0
    par.tmax = stable_time;
else
    par.tmax = -stable_time;
end
par.t0 = par.tmax/1e6;
par.tpoints = pnts;

disp(['Performing initial transient with tmax = ', num2str(par.tmax)]);
sol = df(sol_Rs, par);

% Switch off precondition voltage
par.Vapp = 0;

% Stable time is enetered as negative then run to stable solution
if stable_time < 0
    all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    
    % loop to check ions have reached stable config- if not accelerate ions by
    % order of mag
    j = 1;
    
    while any(all_stable) == 0
        disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10^j)]);
        
        par.tmax = par.tmax*10;
        par.t0 = par.tmax/1e6;
        
        sol = df(sol_Rs, par);
        
        all_stable = verifyStabilization(sol.u, sol.t, 0.7);
        
    end
end

sol_ill = sol;
sol_ill.par.mobseti = par_origin.mobseti;

%dfplot.ELx(sol_ill);
%dfplot.Voct(sol_ill);
disp('LIGHTON complete')

end