function steadystate_struct = stabilize(struct)
% STABILIZE - Simulates the solution increasing the maximum time until when steady state is reached
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
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%% Start code
% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% May 2018; Last revision: May 2018

%------------- BEGIN CODE --------------

% shortcut
par = struct.par;

% set tpoints, just a few ones are necessary here
par.tpoints = 10;
par.tmesh_type = 2; % log spaced time mesh

%% estimate a good tmax
% a tmax too short would make the solution look stable even
% if it's not; too large and the simulation could fail
min_tmax_ions = 10;
min_tmax_freecharges = 1e-3;

% if both mobilities are set
if max(par.mu_c) && par.mu_n(1)
    par.tmax = min([min_tmax_ions, par.tmax*1e4, 2^(-log10(max(par.mu_c))) / 10 + 2^(-log10(par.mu_n(1)))]);
    min_tmax = min_tmax_ions;
% if ionic mobility is zero but free charges mobility is set
elseif par.mu_n(1)
    par.tmax = min([min_tmax_freecharges, par.tmax*1e4, 2^(-log10(par.mu_n(1)))]);
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
warning('off', 'Driftfusion:verifyStabilization');

i=0;
while forceStabilization || ~verifyStabilization(steadystate_struct.u, steadystate_struct.t, 1e-3) % check stability
    i = i + 1;
    if i > 10
        warning('Driftfusion:stabilize', [mfilename ' - not stabile after 10 attempts, giving up.']);
        break
    end

    disp([mfilename ' - Stabilizing ' inputname(1) ' over ' num2str(par.tmax) ' s with an applied voltage of ' num2str(par.Vapp) ' V']);
    % every cycle starts from the last timepoint of the previous cycle
    steadystate_struct = df(steadystate_struct, par);
    if size(steadystate_struct.u, 1) ~= steadystate_struct.par.tpoints % simulation failed
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

% re-enable the warnings
warning('on', 'Driftfusion:verifyStabilization');

%------------- END OF CODE --------------
