function [sol_eq, sol_i_eq, ssol_eq, ssol_i_eq, sol_i_eq_SR] = equilibrate
%EQUILIBRATE Uses analytical initial conditions and runs to equilibrium and steady state
% Takes the parameters from pinParams.m file and tries
% to obtain an equilibrium solution (as if the device has been left for
% a long period of time). This solution can then be used as accurate
% initial conditions for other simulations, e.g. a JV scan.
% Note that tmax is consistently adjusted to appropriate values for to
% ensure there are numerous mesh points where large gradients in the time
% dimension are present.
%
% Syntax:  [sol_eq, sol_eq_SR, sol_i_eq, sol_i_eq_SR, ssol_eq, ssol_eq_SR, ssol_i_eq, ssol_i_eq_SR, sol_light, sol_light_SR, sol_i_light, sol_i_light_SR, ssol_light, ssol_light_SR, ssol_i_light, ssol_i_light_SR] = EQUILIBRATE()
%
% Inputs:
%
% Outputs:
%   sol_eq - short circuit, dark, no mobile ionic defects, no SRH
%   sol_eq_SR - short circuit, dark, no mobile ionic defects, with SRH
%   sol_i_eq - short circuit, dark, mobile ionic defects, no SRH
%   sol_i_eq_SR - short circuit, dark, mobile ionic defects, with SRH
%   ssol_eq - open circuit, dark, no mobile ionic defects, no SRH
%   ssol_eq_SR - open circuit, dark, no mobile ionic defects, with SRH
%   ssol_i_eq - open circuit, dark, mobile ionic defects, no SRH
%   ssol_i_eq_SR - open circuit, dark, mobile ionic defects, with SRH
%   sol_light - short circuit, 1 sun, no mobile ionic defects, no SRH
%   sol_light_SR - short circuit, 1 sun, no mobile ionic defects, with SRH
%   sol_i_light - short circuit, 1 sun, mobile ionic defects, no SRH
%   sol_i_light_SR - short circuit, 1 sun, mobile ionic defects, with SRH
%   ssol_light - open circuit, 1 sun, no mobile ionic defects, no SRH
%   ssol_light_SR - open circuit, 1 sun, no mobile ionic defects, with SRH
%   ssol_i_light - open circuit, 1 sun, mobile ionic defects, no SRH
%   ssol_i_light_SR - open circuit, 1 sun, mobile ionic defects, with SRH
%
% Example:
%   equilibrate()
%     generate stabilized solutions and save them to base workspace
%
% Other m-files required: pindrift, pinParams, mobsetfun
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift, paramsStruct.

% Author: Phil Calado, Ph.D.
% Imperial College London
% Research Group Prof. Jenny Nelson
% email address: p.calado13@imperial.ac.uk
% Contributors: Ilario Gelmetti, Ph.D. student
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% 2015; Last revision: January 2018

%------------- BEGIN CODE --------------

tic;    % Start stopwatch

%% Initial arguments
% Setting sol.sol = 0 enables a parameters structure to be read into
% pindrift but indicates that the initial conditions should be the
% analytical solutions
sol.sol = 0;    

p = pinParams;

% Store initial parameters
original_p = p;

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
disp('Initial solution, zero mobility')
sol = pindrift(sol, p);
disp('Complete')

p.figson = 1;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;

%% Mobility with mobility switched on

% switch on electron and hole mobility
p.mue_i = original_p.mue_i; % electron mobility in intrinsic
p.muh_i = original_p.muh_i; % hole mobility in intrinsic
p.mue_p = original_p.mue_p; % electron mobility in p-type
p.muh_p = original_p.muh_p; % hole mobility in n-type
p.mue_n = original_p.mue_n; % electron mobility in p-type
p.muh_n = original_p.muh_n; % hole mobility in n-type

disp('Solution with mobility switched on')
sol = pindrift(sol, p);

p.Ana = 1;
p.calcJ = 2;
p.tmax = 1e-2;
p.t0 = p.tmax/1e10;

sol_eq = pindrift(sol, p);
disp('Complete')

%% Set up solution for open circuit
disp('Switching boundary conditions to zero flux')
%p.Ana = 0;
p.BC = 0;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;

sol = pindrift(sol_eq, p);
disp('Complete')

%% Symmetricise the solution
disp('Symmetricise solution for open circuit')
symsol = symmetricize(sol);
disp('Complete')

%% Equilibrium solution with mirrored cell and OC boundary conditions, mobility zero
disp('Initial equilibrium open circuit solution')
p.BC = 1;
p.OC = 1;
p.calcJ = 0;
p = mobsetfun(0, 0, p);

ssol = pindrift(symsol, p);
disp('Complete')

%% OC with mobility switched on
disp('Open circuit solution with mobility switched on')
p.tmax = 1e-6;
p.t0 = p.tmax/1e3;

% switch on electron and hole mobility
p.mue_i = original_p.mue_i; % electron mobility in intrinsic
p.muh_i = original_p.muh_i; % hole mobility in intrinsic
p.mue_p = original_p.mue_p; % electron mobility in p-type
p.muh_p = original_p.muh_p; % hole mobility in n-type
p.mue_n = original_p.mue_n; % electron mobility in p-type
p.muh_n = original_p.muh_n; % hole mobility in n-type

ssol = pindrift(ssol, p);

% Longer time step to ensure equilibrium has been reached
p.tmax = 1e-2;
p.t0 = p.tmax/1e3;

ssol_eq = pindrift(ssol, p);
disp('Complete')

%% Equilibrium solutions with ion mobility switched on
%% Closed circuit conditions
disp('Closed circuit equilibrium with ions')

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
disp('Complete')

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
disp('Complete')

% Switch off SR
p.taun_etl = 1e6;
p.taup_etl = 1e6;
p.taun_htl = 1e6;
p.taup_htl = 1e6; 

%% Symmetricise closed circuit condition
disp('Symmetricise equilibriumion solution')
symsol = symmetricize(sol_i_eq);
disp('Complete')

p.OC = 1;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p = mobsetfun(0, 0, p);

%% OC condition with ions at equilbirium
disp('Open circuit equilibrium with ions')
ssol = pindrift(symsol, p);

p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p.mui = 0;

% switch on electron and hole mobility
p.mue_i = original_p.mue_i; % electron mobility in intrinsic
p.muh_i = original_p.muh_i; % hole mobility in intrinsic
p.mue_p = original_p.mue_p; % electron mobility in p-type
p.muh_p = original_p.muh_p; % hole mobility in n-type
p.mue_n = original_p.mue_n; % electron mobility in p-type
p.muh_n = original_p.muh_n; % hole mobility in n-type

ssol = pindrift(ssol, p);

% Switch on ion mobility to ensure equilibrium has been reached
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p.mui = original_p.mui; % this requires mui to be set in the original params

ssol = pindrift(ssol, p);

p.tmax = 1e-2;
p.t0 = p.tmax/1e3;

ssol_i_eq = pindrift(ssol, p);
disp('Complete')

disp('EQUILIBRATION COMPLETE')
toc

end
