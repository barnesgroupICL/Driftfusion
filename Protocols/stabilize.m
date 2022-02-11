function steadystate_struct = stabilize(struct)
% STABILIZE - Simulates the solution increasing the maximum time until when steady state is reached
%
% Syntax:  steadystate_struct = stabilize(struct)
%
% Inputs:
%   STRUCT - a solution struct as created by DF.
%
% Outputs:
%   STEADYSTATE_STRUCT - a solution struct that reached its steady state
%
% Example:
%   soleq_ion_stable = stabilize(soleq.ion);
%     checks if a solution reached steady state and replaces it with its steady state condition
%
% Other m-files required: df, verifyStabilization
% Subfunctions: none
% MAT-files required: none
%
% See also df, verifyStabilization.
%
%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Start code
%------------- BEGIN CODE --------------

% initialize the output as the input, if already stable it will stay like
% this
steadystate_struct = struct;

% shortcut
par_old = struct.par;
par_new = par_old;

% set tpoints, just a few ones are necessary here
par_new.tpoints = 10;
par_new.tmesh_type = 2; % log spaced time mesh

%% estimate a good tmax
% a tmax too short would make the solution look stable even
% if it's not; too large and the simulation could fail

forceStabilization = false;
% checks if the input structure reached the final time point
if size(steadystate_struct.u, 1) ~= par_old.tpoints
    % in case the provided solution did not reach the final time point,
    % force stabilization
    forceStabilization = true;
end
% check if the prefious step was run over sufficient time. In case it was
% not, the solution could seem stable to verifyStabilization just because
% the solution did not evolve much in a time that was very short

min_tmax_freecharges = max((par_old.d.^2)./(2*par_old.kB*par_old.T*min(par_old.mu_n, par_old.mu_p)));

% if both mobilities are set
if par_old.mobseti && par_old.N_ionic_species > 0
    t_c_diff = (par_old.d.^2)./(2*par_old.kB*par_old.T*par_old.mu_c);
    t_c_diff = max(t_c_diff(isfinite(t_c_diff)));
    if par_old.N_ionic_species == 2
        t_a_diff = (par_old.d.^2)./(2*par_old.kB*par_old.T*par_old.mu_a);
        t_a_diff = max(t_a_diff(isfinite(t_a_diff)));
        min_tmax = max(t_c_diff, t_a_diff);
    else
        min_tmax = t_c_diff;
    end
    par_new.tmax = min([min_tmax, min_tmax_freecharges, par_old.tmax*1e4]);
% if ionic mobility is disabled just consider the free charges
% characteristic time
else
    par_new.tmax = min(min_tmax_freecharges, par_old.tmax*1e4);
    min_tmax = min_tmax_freecharges;
end

if par_old.tmax < min_tmax
    forceStabilization = true;
end

par_new.t0 = par_new.tmax / 1e8;

%% equilibrate until steady state

% the warnings are not needed here
warning('off', 'Driftfusion:verifyStabilization');

i=0;
while forceStabilization || ~verifyStabilization(steadystate_struct.u, steadystate_struct.t, 1e-3) % check stability
    i = i + 1;
    if i > 20
        warning('Driftfusion:stabilize', [mfilename ' - not stabile after 20 attempts, giving up.']);
        break
    end

    disp([mfilename ' - Stabilizing ' inputname(1) ' over ' num2str(par_new.tmax) ' s with an applied voltage of ' num2str(getVend(struct)) ' V']);
    % every cycle starts from the last timepoint of the previous cycle
    steadystate_struct = df(steadystate_struct, par_new);
    if size(steadystate_struct.u, 1) ~= steadystate_struct.par.tpoints % simulation failed
        % if the stabilization breaks (does not reach the final time point), reduce the tmax
        par_new.tmax = par_new.tmax / 10;
        par_new.t0 = par_new.tmax / 1e8;
        forceStabilization = true;
    elseif par_new.tmax < min_tmax % did not run yet up to minimum time
        % if the simulation time was not enough, force next step
        par_new.tmax = par_new.tmax * 5;
        par_new.t0 = par_new.tmax / 1e8;
        forceStabilization = true;
    else % normal run
        % each round, increase the simulation time by 5 times
        % (increasing it more quickly can cause some simulation to fail)
        par_new.tmax = par_new.tmax * 5;
        par_new.t0 = par_new.tmax / 1e8;
        forceStabilization = false;
    end

end

% re-enable the warnings
warning('on', 'Driftfusion:verifyStabilization');

%------------- END OF CODE --------------
