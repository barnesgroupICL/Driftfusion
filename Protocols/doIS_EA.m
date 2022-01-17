function struct_IS = doIS_EA(struct_Int, deltaV, freq, periods, tpoints_per_period, stability_timefraction, RelTol)
%DOIS_EA - Do a single Impedance Spectroscopy (IS) simulation
% Starting from the last time point of the provided solution it applies a
% time-varying voltage.
% This voltage profile is the background voltage bias from the provided
% solution summed to a sinusoidally oscillating voltage.
%
% Syntax:  struct_IS = doIS_EA(struct_Int, deltaV, freq, periods, tpoints_per_period, stability_timefraction, RelTol)
%
% Inputs:
%   STRUCT_INT - a single struct as created by DF.
%   DELTAV - voltage oscillation amplitude in volts
%   FREQ - voltage oscillation frequency
%   PERIODS - number of periods to be simulated
%   TPOINTS_PER_PERIOD - how many time points should the solution have for
%     each period, suggested to use a multiple of 4
%   STABILITY_TIMEFRACTION - numeric between 0 and 1. The smaller the
%     stricter the check. When 1: disable the stability check and simply
%     return the first simulation result. When zero or strictly included
%     between 0 and 1: check if the oscillating solution reached a
%     (oscillating) stabilization. In this case, the number rounded to
%     match a whole number of periods and is used as the time_fraction
%     argument in verifyStabilization. Suggested values: 0.5 for normal
%     cases, 0 for cases in which the stabilization is very important,
%     1 for cases in which the first solution is needed.
%   RELTOL - starting relative tolerance of the PDEPE solver, for example a
%     value of 1e-8 sets for a very precise simulation
%
% Outputs:
%   STRUCT_IS - a struct with a solution being perturbed by an
%     oscillating voltage
%
% Example:
%   i_sr_is_100mHz_100mV = doIS_EA(soleq.ion, 0.1, 1e-2, 20, 40, 0.5, 1e-8)
%     simulate an oscillating voltage at 0.1 Hz and 100 mV of half peak to peak voltage amplitude,
%     20 periods and 40 time points per period, taking care of reaching a stable solution,
%     and using a starting relative tolerance of 1e-8, starting from a dark
%     solution at zero applied voltage
%     calculate also the electronic current needed by impedance simulations
%   i_sr_sc_1S_is_100mHz_100mV = doIS_EA(changeLight(soleq.ion, 1, 0), 0.1, 1e-2, 20, 40, 0.5, 1e-8)
%     as above, but with 1 sun illumination starting from short circuit
%     conditions
%   i_sr_oc_1S_is_100mHz_100mV = doIS_EA(findOptimVoc(changeLight(soleq.ion, 1, 0)), 0.1, 1e-2, 20, 40, 0.5, 1e-8)
%     as above, but at starting from open circuit conditions
%
% Other m-files required: df, verifyStabilization, dfana
% Subfunctions: none
% MAT-files required: none
%
% See also df, verifyStabilization, IS_script, EA_script.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

% suggested to use linear time mesh (1), not exponential as usual (2),
% 'cause here we're not interested in a decay, but in a oscillation going on during the whole simulated timespan
logtime = false;

tmax = periods / freq;
% it's good to simulate exactly an integer number of periods, so that the
% last point can be used as the starting point for another simulation,
% useful for reaching stabilization
tpoints = 1 + tpoints_per_period * periods;

% accuracy
struct_Int.par.RelTol = RelTol;

%% define applied voltage and fitting functions
struct_Int.par.JV = 2; % mode for arbitrary Vapp functions

% get Efn and Efp
Vapp_arr = dfana.calcVapp(struct_Int);

% take current voltage, as defined in dfana
Vstart = Vapp_arr(end);

% applied voltage profile, check out the definition of the functions in
% Core/fun_gen.m
Vapp_func = 'sin';
Vapp_coeff = [Vstart, deltaV, freq, 0];

%% first run
disp([mfilename ' - Int: ' num2str(struct_Int.par.int1) '; Vdc: ' num2str(Vapp_arr(end)) ' V; Freq: ' num2str(freq) ' Hz']);
struct_IS = VappFunction(struct_Int, Vapp_func, Vapp_coeff, tmax, tpoints, logtime);

%% run until oscillating stabilization
% verify that the final time point has the same solution as the last point
% of an oscillation in the middle of the simulated timespan, so verify that "periods" was enough
% if just the first solution is requested setting reach_stability to
% false, break after the first cycle
i = 0;
warning('off', 'Driftfusion:verifyStabilization'); % verifyStabilization warnings are substituted by others here
while ~verifyStabilization(struct_IS.u, struct_IS.t, floor(stability_timefraction * periods) / periods)
    disp([mfilename ' - Freq: ' num2str(freq) ' Hz solution was not stabilized, trying again'])
    struct_IS = VappFunction(struct_IS, Vapp_func, Vapp_coeff, tmax, tpoints, logtime);
    i = i+1;
    % in case 10 repetitions were not enough, stop
    if i > 10
        warning('Driftfusion:IS_EA',...
            'IS_EA seems that the solution did not reach complete stabilization after %s repetitions on %s periods',...
            num2str(i), num2str(periods))
        break
    end
end
warning('on', 'Driftfusion:verifyStabilization'); % re-enable warnings

%------------- END OF CODE --------------
