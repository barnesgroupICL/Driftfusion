function sol_int = changeLight(sol, newInt, tmax, lightSource)
% CHANGELIGHT - Stabilize solutions at a new light intensity
%
% Syntax:  sol_int = changeLight(sol, newInt, tmax, lightSource)
%
% Inputs:
%   SOL - a solution sol as created by DF.
%   NEWINT - float, the requested light intensity, both increasing,
%     decreasing and setting to zero are supported (but for a solution at
%     zero illumination is more meaningful to take the solution directly
%     from equilibrate)
%   TMAX - float, the time over which simulate the illumination change,
%     can be zero for an automatic guess
%   LIGHTSOURCE - optional integer, illumination source to be varied:
%     1 = first illumination source (from sol.par.light_source1),
%     2 = second illumination source (from sol.par.light_source2)
%
% Outputs:
%   sol_int - a solution sol at NEWINT light intensity
%
% Example:
%   sol_ion_1S = changeLight(soleq.ion, 1, 0)
%     take the solution soleq.ion and stabilize to a new light intensity
%     of 1 sun, guess a good tmax time, the light source is not specified,
%     so the first one is used (usually, a solar spectra)
%   sol_ion_5lazor = changeLight(soleq.ion, 5, 0, 2)
%     as above, but set the intensity of the secondary illumination source
%     (usually, a monochromatic illumination)
%

% Other m-files required: df, stabilize
% Subfunctions: none
% MAT-files required: none
%
% See also genIntStructs, df.
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
% October 2017; Last revision: January 2021

%------------- BEGIN CODE --------------

par = sol.par;
par.tmesh_type = 2;
par.t0 = 1e-10;
par.tpoints = 30;

% set an initial time for the simulation
if tmax
    par.tmax = tmax;
else % if tmax was zero, estimate a good one
    if par.mobseti
        tmax_temp = min(1, 2^(-log10(par.mu_c(par.active_layer(1)))) / 10 + 2^(-log10(par.mu_n(par.active_layer(1)))));
    else
        tmax_temp = min(1e-3, 2^(-log10(par.mu_n(par.active_layer(1)))));
    end
end

if ~exist('lightSource','var')
    lightSource = 1;
end

oldInt = getInt(par, lightSource);

%% change light in steps
% not needed to reach a good stabilized solution in
% each step, so stabilization is not verified here

lastStepSuccessful = false;
steps = 1;

while steps < 7 && ~lastStepSuccessful
    sol_int = sol;
    % list of valid illumination values to use
    % sum a small number to avoid the problems with zero illumination
    tiny = max(oldInt,newInt)/10000;
    intArray = logspace(log10(oldInt+tiny), log10(newInt+tiny), steps+1);
    % remove oldInt from the list
    intArray = intArray(2:end)-tiny;
    for i = intArray
        prevInt = getInt(sol_int.par, lightSource);
        switch lightSource
            case 1
                par.g1_fun_type = 'sweepAndStill';
                % COEFF = [A_start, A_end, sweep_duration]
                par.g1_fun_arg = [prevInt, i, par.tmax/2e3];
            case 2
                par.g2_fun_type = 'sweepAndStill';
                % COEFF = [A_start, A_end, sweep_duration]
                par.g2_fun_arg = [prevInt, i, par.tmax/2e3];
        end
        disp([mfilename ' - Go from light intensity ' num2str(prevInt) ' to ' num2str(i) ' over ' num2str(par.tmax) ' s'])

        sol_int = df(sol_int, par);

        if size(sol_int.u,1) == par.tpoints
            lastStepSuccessful = true;
        else
            lastStepSuccessful = false;
            steps = steps + 1;
            break
        end
    end
end

if ~lastStepSuccessful
    error('Driftfusion:changeLight', 'THE SOLUTION DID NOT REACH COMPLETION')
end

%at some point in the future par.int1 should be eliminated in favour of g1_fun_arg(1)
switch lightSource
    case 1
        sol_int.par.g1_fun_type = 'constant';
        sol_int.par.g1_fun_arg = newInt;
        if isprop(sol_int.par, 'int1')
            sol_int.par.int1 = newInt;
        end
    case 2
        sol_int.par.g2_fun_type = 'constant';
        sol_int.par.g2_fun_arg = newInt;
        if isprop(sol_int.par, 'int2')
            sol_int.par.int2 = newInt;
        end
end

%% stabilize
sol_int = stabilize(sol_int); % go to steady state

function int = getInt(par, lightSource)
%at some point in the future par.int1 should be eliminated in favour of g1_fun_arg(1)
switch lightSource
    case 1
        if isprop(par, 'int1')
            int = par.int1;
        else
            int = par.g1_fun_arg(1);
        end
    case 2
        if isprop(par, 'int2')
            int = par.int2;
        else
            int = par.g2_fun_arg(1);
        end
end


%------------- END OF CODE --------------
