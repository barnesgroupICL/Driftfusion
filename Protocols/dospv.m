function spvsol = dospv(sol_ini, Int, mobseti, tpoints, tmax, Rs, stabilise)
% A function for switching the light on and off with a series
% resistance Rs

%% Input argumanets
% SOL_INI is the initial solution (usually equilibrium)
% INT = Light intensity in Suns
% MOBSETI = 1 or 0 - mobility prefactor can be used to switch off ion
% mobility
% TPOINTS = number of points in the time array- 200 is good
% TMAX = Initial guess for the end time. 
% RS = Series resistance. 1e6 Ohms will give you approximately open circuit
% STABILISE = Choose whether to reach a steady-state or not. If STABILISE = 1 then the code
% will loop until stable

disp('Starting SPV')

% Get initial parameters
par = sol_ini.par;

%% Set parameters for initial dark solution with Rs on
par.mobseti = 0;
par.Int = 0;
par.tmesh_type = 1;
par.tpoints = 100;
par.tmax = 1e-3;
par.t0 = 0;
par.Rs = Rs;
par.Rs_initial = 1;

%% Get dark solution with the series resistance switched on- this is
% necessary because the current at equilibrium does not reach zero. This
% step applies a small voltage simulating the series resistance
disp(['initial dark solution with Rs = ', num2str(Rs), ' Ohms switched on'])

sol_Rs = df(sol_ini, par);

disp('Complete')
%% Parameters for initial illuminated solution
par.Rs_initial = 0;
par.int1 = Int;
par.tmesh_type = 2;
par.tpoints = tpoints;
par.tmax = tmax;
par.t0 = par.tmax/1e6;
par.mobseti = mobseti;

%% Get initial illuminated solution
disp(['SPV initial solution, tmax = ', num2str(par.tmax)])
sol = df(sol_Rs, par);
disp('Complete')

%% Run stabilisation if selected
if stabilise
    all_stable = verifyStabilization(sol.u, sol.t, 0.9);
    
    j = 1;
    while any(all_stable) == 0 && j<6
        
        par.tmax = par.tmax*10;
        par.t0 = par.tmax/1e6;
        disp(['increasing SPV run time, tmax = ', num2str(par.tmax)]);
        
        sol = df(sol_ini, par);
        
        all_stable = verifyStabilization(sol.u, sol.t, 0.9);
        j = j+1;
    end
end

%% store the end solution
sol_ill = sol;

%% Clear parameters
spvsol.ill = sol_ill;
spvsol.ill.par.mobseti = 1;
spvsol.ill.par.figson = 1;

%% Dark step
par.int1 = 0;

disp('SPV dark step')
sol_dk = df(sol_ill, par);

spvsol.dk = sol_dk;
spvsol.dk.par.mobseti = 1;
spvsol.dk.par.figson = 1;
spvsol.Int = Int;

%spvsol.dat = spvana5(spvsol);

disp('SPV complete')

end



