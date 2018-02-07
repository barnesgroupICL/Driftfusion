function [sol_eq, sol_i_eq, ssol_eq, ssol_i_eq, sol_i_eq_SR] = equilibrate
% Uses analytical initial conditions and runs to equilibrium
% Note that tmax is consistently adjusted to appropriate values for to
% ensure there are numerous mesh points where large gradients in the time
% dimension are present

tic;    % Start stopwatch

%% Initial arguments
% Setting sol.sol = 0 enables a parameters structure to be read into
% pindrift but indicates that the initial conditions should be the
% analytical solutions
sol.sol = 0;    

p = pinParams;

% Store initial mobility of intrinsic layer- note all mobilities will be
% set to this value during the equilibration procedure.
mue_i = p.mue_i;

%% Start with low recombination coefficients
p.klin = 0;
p.klincon = 0;
p.taun_etl = 1e6;       % [s] SRH time constant for electrons
p.taup_etl = 1e6;      % [s] SRH time constant for holes
p.taun_htl = 1e6;       %%%% USE a high value of (e.g.) 1 to switch off
p.taup_htl = 1e6;

% Raditative recombination could also be set to low values initially if required. 
% p.krad = 1e-20;
% p.kradetl = 1e-20;
% p.kradhtl = 1e-20;

%% General initial parameters
p.tmesh_type = 2;
p.tpoints = 200;

p.Ana = 0;
p.JV = 0;
p.Vapp = 0;
p.Int = 0;
p.pulseon = 0; 
p.OC = 0;
p.BC = 1;
p.tmesh_type = 2;
p.tmax = 1e-9;
p.t0 = p.tmax/1e4;

%% Mobsetfun is used to easily set mobilities
p = mobsetfun(0, 0, p);

%% Initial solution with zero mobility
disp('initial solution, zero mobility')
sol = pindrift(sol, p);
disp('complete')

p = mobsetfun(mue_i, 0, p);
p.figson = 1;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;

%% Mobility with mobility switched on
disp('solution with mobility switched on')
sol = pindrift(sol, p);

p.Ana = 1;
p.calcJ = 2;
p.tmax = 1e-2;
p.t0 = p.tmax/1e10;

sol_eq = pindrift(sol, p);
disp('complete')

%% Set up solution for open circuit
disp('switch boundary conditions to zero flux')
%p.Ana = 0;
p.BC = 0;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;

sol = pindrift(sol_eq, p);
disp('complete')

%% Symmetricise the solution
disp('symmetricise solution for open circuit')
symsol = symmetricize(sol);
disp('complete')

%% Equilibrium solution with mirrored cell and OC boundary conditions, mobility zero
disp('initial equilibrium open circuit solution')
p.BC = 1;
p.OC = 1;
p.calcJ = 0;
p = mobsetfun(0, 0, p);

ssol = pindrift(symsol, p);
disp('complete')

%% OC with mobility switched on
disp('open circuit solution with mobility switched on')
p.tmax = 1e-6;
p.t0 = p.tmax/1e3;
p = mobsetfun(mue_i, 0, p);

ssol = pindrift(ssol, p);

% Longer time step to ensure equilibrium has been reached
p.tmax = 1e-2;
p.t0 = p.tmax/1e3;
p = mobsetfun(mue_i, 0, p);

ssol_eq = pindrift(ssol, p);
disp('complete')

%% Equilibrium solutions with ion mobility switched on
%% Closed circuit conditions
disp('closed circuit equilibrium with ions')

p.OC = 0;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p.mui = 1e-6;           % Ions are accelerated to reach equilibrium

sol = pindrift(sol_eq, p);

% Much longer second step to ensure that ions have migrated
p.calcJ = 2;
p.tmax = 1e-2;
p.t0 = p.tmax/1e3;

sol_i_eq = pindrift(sol, p);
disp('complete')

%% Ion equilibrium with surface recombination
disp('Switching on surface recombination')
p.taun_etl = 1e-11;
p.taup_etl = 1e-11;
p.taun_htl = 1e-11;
p.taup_htl = 1e-11; 

p.calcJ = 0;
p.tmax = 1e-6;
p.t0 = p.tmax/1e3;

sol_i_eq_SR = pindrift(sol_i_eq, p);
disp('complete')

% Switch off SR
p.taun_etl = 1e6;
p.taup_etl = 1e6;
p.taun_htl = 1e6;
p.taup_htl = 1e6; 

%% Symmetricise closed circuit condition
disp('symmetricise equilibriumion solution')
symsol = symmetricize(sol_i_eq);
disp('complete')

p.OC = 1;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p = mobsetfun(0, 0, p);

%% OC condition with ions at equilbirium
disp('open circuit equilibrium with ions')
ssol = pindrift(symsol, p);

p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p.mui = 0;
p = mobsetfun(mue_i, 0, p);

ssol = pindrift(ssol, p);

% Switch on ion mobility to ensure equilibrium has been reached
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p = mobsetfun(mue_i, 1e-6, p);

ssol = pindrift(ssol, p);

p.tmax = 1e-2;
p.t0 = p.tmax/1e3;

ssol_i_eq = pindrift(ssol, p);
disp('complete')

disp('Equilibration complete')
toc

end
