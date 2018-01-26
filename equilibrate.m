function [sol_eq, sol_light, ssol_eq, ssol_light, sol_i_eq, sol_i_light, ssol_i_eq, ssol_i_light] = equilibrate(name)
%EQUILIBRATE Generate accurate initial conditions at equilibrium
%Takes the parameters from pinParams.m file and tries
%to obtain an equilibrium solution (as if the device has been left for
%a long period of time). This solution can then be used as accurate
%initial conditions for other simulations, e.g. a JV scan.
%
% Syntax:  [sol_eq, sol_light, ssol_eq, ssol_light, sol_i_eq, sol_i_light, ssol_i_eq, ssol_i_light] = EQUILIBRATE(NAME)
%
% Inputs:
%   NAME - if is provided, saves the solution as [NAME *sol*] in the
%   volatile base workspace.
%
% Outputs:
%   sol_eq - short circuit in dark without mobile ionic defects
%   sol_light - short circuit at 1 sun conditions without mobile ionic defects
%   ssol_eq - open circuit in dark without mobile ionic defects
%   ssol_light - open circuit at 1 sun conditions without mobile ionic defects
%   sol_i_eq - short circuit in dark with mobile ionic defects
%   sol_i_light - short circuit at 1 sun conditions with mobile ionic defects
%   ssol_i_eq - open circuit in dark without mobile ionic defects
%   ssol_i_light - open circuit at 1 sun conditions without mobile ionic defects
%
% Example:
%   Equilibrate()
%     generate stabilized solutions and save them to base workspace
%   Equilibrate('new_')
%     as above but prepending a string to the saved struct name
%
% Other m-files required: pindrift, pinParams, mobsetfun
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift, paramsStruct.

% Author: Phil Calado, Ph.D., Ilario Gelmetti, Ph.D. student
% Imperial College London
% Research Group Prof. Jenny Nelson
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: p.calado13@imperial.ac.uk, iochesonome@gmail.com
% Supervised by: Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% 2015; Last revision: January 2018

%------------- BEGIN CODE --------------

if ~nargin
    name = '';
end

% initial empty solution makes pindrift to generate an analytical one
sol.sol = 0;

% load parameters from pinParams.m file
p = pinParams;
original_p = p;

p.tmesh_type = 2; % logarithmic time mesh
p.tpoints = 30; % few points, just the last one is usually used
p.JV = 0; % no voltage scan
p.Vapp = 0; % no applied voltage
p.Int = 0; % dark
p.pulseon = 0; % no transient illumination pulse
p.OC = 0; % short circuit
p.calcJ = 0; % no need for current calculation

%% SC condition without ions in dark

% freeze electrons and holes and ions for an initial solution
p = mobsetfun(0, 0, p);

p.tmax = 1e-9;
p.t0 = p.tmax/1e3;

sol_eq = pindrift(sol, p);

% switch on electron and hole mobility
p.mue_i = original_p.mue_i; % electron mobility in intrinsic
p.muh_i = original_p.muh_i; % hole mobility in intrinsic
p.mue_p = original_p.mue_p;
p.muh_p = original_p.muh_p;
p.mue_n = original_p.mue_n;
p.muh_n = original_p.muh_n;

p.tmax = 1;
p.t0 = p.tmax/1e9;

sol_eq = pindrift(sol_eq, p);
assignin('base', [name 'sol_eq'], sol_eq);

%% SC condition without ions at 1 sun conditions

p = sol_eq.params;
p.Int = 1;

sol_light = pindrift(sol_eq, p);
assignin('base', [name 'sol_light'], sol_light);

%% OC condition without ions at equilibrium

symsol = symmetricize(sol_eq);

p = sol_eq.params;
p.OC = 1;

ssol_eq = pindrift(symsol, p);
assignin('base', [name 'ssol_eq'], ssol_eq);

%% OC condition without ions at 1 sun conditions

p = ssol_eq.params;
p.Int = 1;

ssol_light = pindrift(ssol_eq, p);
assignin('base', [name 'ssol_light'], ssol_light);

%% SC condition with ions in dark

p = sol_eq.params;

%Switch on ion mobility
p.mui = original_p.mui; % ionic mobility

sol_i_eq = pindrift(sol_eq, p);
assignin('base', [name 'sol_i_eq'], sol_i_eq);

%% SC condition with ions at 1 sun conditions

p = sol_i_eq.params;
p.Int = 1;

sol_i_light = pindrift(sol_i_eq, p);
assignin('base', [name 'sol_i_light'], sol_i_light);

%% OC condition with ions in dark

symsol = symmetricize(sol_i_eq);

p = sol_i_eq.params;
p.OC = 1;

warning('off', 'pindrift:verifyStabilization'); % the following line is expected to fail on stabilization
ssol_i_eq = pindrift(symsol, p); % first run is not enough for stabilization

% ssol_i_eq takes long time to equilibrate
p.tmax = 2 ^ (-log10(p.mui)) + 1e4 * 2 ^ (-log10(p.mue_i));
p.t0 = p.tmax / 1e9;

ssol_i_eq = pindrift(ssol_i_eq, p);
warning('on', 'pindrift:verifyStabilization'); % reenable the warnings

ssol_i_eq = pindrift(ssol_i_eq, p);
assignin('base', [name 'ssol_i_eq'], ssol_i_eq);

%% OC condition with ions at 1 sun conditions

p = ssol_i_eq.params;
p.Int = 1;

p.tmax = 10;
p.t0 = p.tmax / 1e9;

warning('off', 'pindrift:verifyStabilization'); % the following line is expected to fail on stabilization
ssol_i_light = pindrift(ssol_i_eq, p);
warning('on', 'pindrift:verifyStabilization'); % reenable the warnings

p.tmax = 2 ^ (-log10(p.mui)) + 1e4 * 2 ^ (-log10(p.mue_i));
ssol_i_light = pindrift(ssol_i_light, p);

assignin('base', [name 'ssol_i_light'], ssol_i_light);

evalin('base', 'clear symsol')
evalin('base', 'clear sol')
evalin('base', 'clear ssol')

%------------- END OF CODE --------------
