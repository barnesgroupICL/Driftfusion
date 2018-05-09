function steadystate_struct = stabilize(struct)
%STABILIZE - Simulates the solution increasing the maximum time until when steady state is reached
%
% Syntax:  steadystate_struct = stabilize(struct)
%
% Inputs:
%   STRUCT - a solution struct as created by PINDRIFT.
%
% Outputs:
%   STEADYSTATE_STRUCT - a solution struct that reached its steady state
%
% Example:
%   ssol_i_1S_SR = stabilize(ssol_i_1S_SR);
%     checks if a solution reached steady state and replaces it with its steady state condition
%
% Other m-files required: pindrift, verifyStabilization
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift, verifyStabilization.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% May 2018; Last revision: May 2018

%------------- BEGIN CODE --------------

% eliminate JV configuration
struct.p.JV = 0;

% shortcut
p = struct.p;

% set tpoints, just a few ones are necessary here
p.tpoints = 5;
p.tmesh_type = 2; % log spaced time mesh

%% estimate a good tmax
% a tmax too short would make the solution look stable even
% if it's not; too large and the simulation could fail

% if both mobilities are set
if p.mui && p.mue_i
    p.tmax = min(5, 2^(-log10(p.mui)) / 10 + 2^(-log10(p.mue_i)));
% if ionic mobility is zero but free charges mobility is set
elseif p.mue_i
    p.tmax = min(1e-3, 2^(-log10(p.mue_i)));
% if no mobility (or just the ionic one) is set
else
    p.tmax = 1e-6;
end

p.t0 = p.tmax / 1e8;

%% obtain the applied voltage Vapp

% in case of open circuit solution, the Vapp is zero
if p.OC
    p.Vapp = 0;
% otherwise use the existing voltage as Vapp
else
    % taken from pinAna, using the last timepoint
    n = struct.sol(end, :, 1); % electrons
    P = struct.sol(end, :, 2); % holes
    V = struct.sol(end, :, 4); % electric potential
    % Calculate energy levels and chemical potential         
    V = V - p.EA;
    % take the right extrema for electrons (boundary for ETL) and left
    % extrema for the holes (boundaryv for HTL)
    Efn_right = real(-V(end) + p.Ei + (p.kB*p.T/p.q) * log(n(end)/p.ni)); % Electron quasi-Fermi level at right boundary
    Efp_left = real(-V(1) + p.Ei - (p.kB*p.T/p.q) * log(P(1)/p.ni)); % Hole quasi-Fermi level at left boundary
    p.Vapp = Efn_right - Efp_left;
end

%% equilibrate until steady state

% initialize the output as the input, if already stable it will stay like
% this
steadystate_struct = struct;

% the warnings are not needed here
warning('off', 'pindrift:verifyStabilization');

while ~verifyStabilization(steadystate_struct.sol, steadystate_struct.t, 1e-4) % check stability
    disp([mfilename ' - Stabilizing ' inputname(1) ' over ' num2str(p.tmax) ' s with an applied voltage of ' num2str(p.Vapp) ' V']);
    % every cycle starts from the last timepoint of the previous cycle
    steadystate_struct = pindrift(steadystate_struct, p);
    % each round, increase the simulation time by 5 times
    % (increasing it more quickly can cause some simulation to fail)
    p.tmax = p.tmax * 5;
end

% re-enable the warnings
warning('on', 'pindrift:verifyStabilization');

%------------- END OF CODE --------------
