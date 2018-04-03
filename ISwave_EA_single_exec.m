function asymstruct_ISwave = ISwave_EA_single_exec(asymstruct_Int, BC, deltaV, freq, periods, tpoints_per_period, reach_stability, EA, RelTol)
%ISWAVE_SINGLE_EXEC - Do a single Impedance Spectroscopy (ISwave) experiment
%
% Syntax:  asymstruct_ISwave = ISwave_EA_single_exec(asymstruct_Int, BC, deltaV, freq, periods, tpoints_per_period, reach_stability, EA, RelTol)
%
% Inputs:
%   ASYMSTRUCT_INT - a single asymmetric struct as created by PINDRIFT.
%   BC - boundary conditions during voltage oscillation
%   DELTAV - voltage oscillation amplitude in volts
%   FREQ - voltage oscillation frequency
%   PERIODS - number of periods to be simulated
%   TPOINTS_PER_PERIOD - how many time points should the solution have for
%     each period, suggested to use a multiple of 4
%   REACH_STABILITY - logical, check if the oscillating solution reached a
%     (oscillating) stabilization, otherwise just use the result of the
%     initial simulation. This can be useful when it's known that the
%     starting solution is not stabilized, for example measuring with an
%     unstabilized ionic profile without frozen_ions
%   EA - logical, set to true when performing ElectroAbsorbance simulation,
%     this will skip the electrical current calculation
%   RELTOL - starting relative tolerance of the PDEPE solver
%
% Outputs:
%   ASYMSTRUCT_ISWAVE - a struct with a solution being perturbated by an
%     oscillating voltage
%
% Example:
%   ISwave_EA_single_exec(asymmetricize(ssol_i_light, 1), 1, 2e-3, 1e6, 20, 40, true, true, 1e-6)
%     simulate an oscillating voltage at 1 MHz and 2 mV of amplitude, with
%     selective contacts, 20 periods and 40 time points per period,
%     calculating the ionic current and using a starting relative tolerance
%     of 1e-4
%
% Other m-files required: pindrift, verifyStabilization
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift, verifyStabilization, ISwave_full_exec.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: March 2018

%------------- BEGIN CODE --------------

p = asymstruct_Int.p;

% suggested to use linear time mesh (1), not exponential as usual (2),
% 'cause here we're not interested in a decay, but in a oscillation going on during the whole simulated timespan
p.tmesh_type = 1;

p.tmax = periods / freq;
p.tpoints = 1 + tpoints_per_period * periods;

if EA
    p.calcJ = 0; % in ElectroAbsorbance simulation we don't need the electrical current value
else
    p.calcJ = 1; % in Impedance Spectroscopy we need the current at the electrodes
    p.Ana = 1;
end

p.BC = BC; % usually BC 1 is ok, 2 can be used
p.RelTol = RelTol; % in case some were defined, for increasing accuracy

%% define applied voltage and fitting functions
p.JV = 2; % new mode for arbitrary Vapp functions

% get Efn and Efp
[~, ~, ~, Efn, Efp, ~] = pinAna(asymstruct_Int);

Vstart = Efn(end, end) - Efp(end, 1);
p.Vapp_params = [Vstart, deltaV, 0, 2 * pi * freq];
p.Vapp_func = @(coeff, t) coeff(1) + coeff(2) * sin(coeff(3) + coeff(4) * t);
p.J_E_func = @(coeff, t) coeff(1) + coeff(2) * sin(coeff(3) + 2 * pi * freq * t);
p.J_E_func_tilted = @(coeff, t, tilting, t_middle) p.J_E_func(coeff, t) + tilting * (t - t_middle);
% double frequency for E squared fitting, useful for Stark spectroscopy
% (ElectroAbsorption EA) simulation
p.E2_func = @(coeff, t) coeff(1) + coeff(2) * sin(coeff(3) + 4 * pi * freq * t);

%% first run
disp([mfilename ' - Int: ' num2str(p.Int) '; Vdc: ' num2str(Efn(end, end) - Efp(end, 1)) ' V; Freq: ' num2str(freq) ' Hz']);
warning('off', 'pindrift:verifyStabilization'); % an oscillating solution is not going to be stable ever
asymstruct_ISwave = pindrift(asymstruct_Int, p);

%% run until stabilization
% verify that the final time point has the same solution as the last point
% of the previous previous oscillation, so verify that "periods" was enough
% if just the first solution is requested setting reach_stability to
% false, break after the first cycle
i = 0;
while ~verifyStabilization(asymstruct_ISwave.sol, asymstruct_ISwave.t, (periods - 2) / periods) && reach_stability
    disp([mfilename ' - solution was not stabilized, trying again'])
    asymstruct_ISwave = pindrift(asymstruct_ISwave, p);
    i = i+1;
    if i > 2
        warning('pindrift:ISwave_EA_single_exec',...
            'ISwave_EA_single_exec seems that the solution did not reach complete stabilization after %s repetitions on %s periods',...
            num2str(i), num2str(periods))
        break
    end
end
warning('on', 'pindrift:verifyStabilization'); % re-enable warnings

%------------- END OF CODE --------------
