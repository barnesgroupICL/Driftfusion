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
% Other m-files required: df, verifyStabilization
% Subfunctions: none
% MAT-files required: none
%
% See also df, verifyStabilization.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% May 2018; Last revision: May 2018

%------------- BEGIN CODE --------------

% eliminate JV configuration
struct.par.JV = 0;

% shortcut
par = struct.par;

% set tpoints, just a few ones are necessary here
par.tpoints = 10;
par.tmesh_type = 2; % log spaced time mesh

% disable Ana, as the ploting done in pinAna results in an error if the simulation doesn't reach the final time point
par.Ana = 0;

%% estimate a good tmax
% a tmax too short would make the solution look stable even
% if it's not; too large and the simulation could fail

min_tmax_ions = 10;
min_tmax_freecharges = 1e-3;

% if both mobilities are set
if max(par.mucat) && par.mue(1)
    par.tmax = min([min_tmax_ions, par.tmax*1e4, 2^(-log10(max(par.mucat))) / 10 + 2^(-log10(par.mue(1)))]);
    min_tmax = min_tmax_ions;
% if ionic mobility is zero but free charges mobility is set
elseif par.mue(1)
    par.tmax = min([min_tmax_freecharges, par.tmax*1e4, 2^(-log10(par.mue(1)))]);
    min_tmax = min_tmax_freecharges;
end

par.t0 = par.tmax / 1e8;

%% equilibrate until steady state

% initialize the output as the input, if already stable it will stay like
% this
steadystate_struct = struct;

forceStabilization = false;
% checks if the input structure reached the final time point
if size(steadystate_struct.u, 1) ~= steadystate_struct.par.tpoints
    % in case the provided solution did not reach the final time point,
    % force stabilization
    forceStabilization = true;
end
% check if the prefious step was run over sufficient time. In case it was
% not, the solution could seem stable to verifyStabilization just because
% the solution did not evolve much in a time that was very short
if struct.par.tmax < par.tmax
    forceStabilization = true;
end

% the warnings are not needed here
warning('off', 'df:verifyStabilization');

while forceStabilization || ~verifyStabilization(steadystate_struct.sol, steadystate_struct.t, 1e-3) % check stability
    disp([mfilename ' - Stabilizing ' inputname(1) ' over ' num2str(par.tmax) ' s with an applied voltage of ' num2str(par.Vapp) ' V']);
    % every cycle starts from the last timepoint of the previous cycle
    steadystate_struct = df(steadystate_struct, par);
    if size(steadystate_struct.sol, 1) ~= steadystate_struct.par.tpoints % simulation failed
        % if the stabilization breaks (does not reach the final time point), reduce the tmax
        par.tmax = par.tmax / 10;
        par.t0 = par.tmax / 1e8;
        forceStabilization = true;
    elseif par.tmax < min_tmax % did not run yet up to minimum time
        % if the simulation time was not enough, force next step
        par.tmax = min(par.tmax * 5, 1e4);
        par.t0 = par.tmax / 1e8;
        forceStabilization = true;
    else % normal run
        % each round, increase the simulation time by 5 times
        % (increasing it more quickly can cause some simulation to fail)
        par.tmax = min(par.tmax * 5, 1e4);
        par.t0 = par.tmax / 1e8;
        forceStabilization = false;
    end

end

% restore original Ana value
steadystate_struct.par.Ana = struct.par.Ana;

% re-enable the warnings
warning('on', 'df:verifyStabilization');

%------------- END OF CODE --------------
