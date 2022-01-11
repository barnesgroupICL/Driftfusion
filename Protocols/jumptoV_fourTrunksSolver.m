function sol_dwell = jumptoV_fourTrunksSolver(sol_ini, Vjump, tdwell, tjump, maxDepth, debug)
%JUMPTOV_FOURTRUNKSSOLVER - an alternative script for jumptoV that tries harder to solve over long dwell times
%
% Syntax:  sol_dwell = jumptoV_fourTrunksSolver(sol_ini, Vjump, tdwell, tjump, maxDepth, debug)
%
% Inputs:
%   SOL_INI - a single struct created by DF.
%   VJUMP - a float, the amplitude of the voltage jump to be performed
%     before the dwell step.
%   TDWELL - a float, the dwell time to evolve the solution at a constant
%     voltage after the jump.
%   TJUMP - optional float, the time used for increasing the voltage.
%   MAXDEPTH - optional integer, the number of iterations the runDfFourTrunks
%     solver could use when the PDEPE solver breaks. An Inf value can be
%     provided for removing the limit.
%   DEBUG - optional logical, set true if runDfFourTrunks should save to
%     the workspace all the partial solutions.
%   
%
% Outputs:
%   sol_dwell - a DF solution structure
%
% Example:
%   soleq_ion_10V_100ks = jumptoV_fourTrunksSolver(soleq.ion, 10, 1e5);
%     jump to 10 V and then evolve the solution for 100'000 seconds. Beware
%     that the output can have an unexpected number of timepoints,
%     depending on how the solving goes.
%
% Other m-files required: df, runDfFourTrunks
% Subfunctions: none 
% MAT-files required: none
%
% See also df, jumptoV, runDfFourTrunks.
%
%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%------------- BEGIN CODE --------------

if tdwell <= 0
    error('please choose a tdwell value > 0 s')
end

%% OUTPUTS
% SOL_RELAX - the output relaxation solution
par = sol_ini.par;

%% Set up sweep
% Characteristic diffusion time
if nargin < 6
    debug = false;
    if nargin < 5
        maxDepth = Inf;
    end
end
t_diff = max((par.d.^2)./(2*par.kB*par.T*min(par.mu_n, par.mu_p)));
t_cationdiff = (par.d.^2)./(2*par.kB*par.T*par.mu_c);
t_cationdiff = max(t_cationdiff(isfinite(t_cationdiff)));

if nargin < 4
    tmax_sweep = t_diff;
else
    tmax_sweep = tjump;
end
par.tmax = tmax_sweep;
par.tmesh_type = 1;
par.t0 = 0;
par.V_fun_type = 'sweep';
par.V_fun_arg(1) = par.Vapp;
par.V_fun_arg(2) = Vjump;
par.V_fun_arg(3) = tmax_sweep;
par.tpoints = 100;

disp([mfilename ' - Jump to ' num2str(Vjump) ' V over ' num2str(tmax_sweep) ' s'])
sol_jump = df(sol_ini, par);
if debug
    assignin('base', "sol_jump", sol_jump)
end

assert(size(sol_jump.u,1) == sol_jump.par.tpoints, [mfilename ' - The jump solution did not get to ' num2str(Vjump) ' V! Try with a different tjump or Vjump.'])

par = sol_jump.par;
disp([mfilename ' - Jump over ' num2str(par.tmax) ' s completed, starting dwell'])

par.V_fun_type = 'constant';
par.V_fun_arg = par.Vapp;

par.tmax = tdwell;
par.t0 = tmax_sweep;
par.tmesh_type = 2;

%sol_dwell = df(sol_jump, par);
%if size(sol_dwell.u,1) ~= par.tpoints
    disp([mfilename ' - Breaking the solution in four trunks'])
    suggested_times = [10*t_diff, 1e3*t_diff, min(1e5*t_diff, t_cationdiff)];
    par.RelTol=1e-4;
    sol_dwell = runDfFourTrunks(sol_jump, par, maxDepth, '', 0, debug, suggested_times);
%end

end

%------------- END OF CODE --------------