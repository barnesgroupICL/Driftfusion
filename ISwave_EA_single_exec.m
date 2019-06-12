function asymstruct_ISwave = ISwave_EA_single_exec(asymstruct_Int, deltaV, freq, periods, tpoints_per_period, reach_stability, EA, RelTol)
%ISWAVE_EA_SINGLE_EXEC - Do a single Impedance Spectroscopy (ISwave) simulation
% Starting from the last time point of the provided solution it applies a
% time-varying voltage.
% This voltage profile is the background voltage bias from the provided
% solution summed to a sinusoidally oscillating voltage.
%
% Syntax:  asymstruct_ISwave = ISwave_EA_single_exec(asymstruct_Int, deltaV, freq, periods, tpoints_per_period, reach_stability, EA, RelTol)
%
% Inputs:
%   ASYMSTRUCT_INT - a single asymmetric struct as created by PINDRIFT.
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
%   RELTOL - starting relative tolerance of the PDEPE solver, for example a
%     value of 1e-8 sets for a very precise simulation
%
% Outputs:
%   ASYMSTRUCT_ISWAVE - a struct with a solution being perturbed by an
%     oscillating voltage
%
% Example:
%   sol_i_eq_SR_is_100mHz_2mV = ISwave_EA_single_exec(sol_i_eq_SR, 2e-3, 1e-2, 20, 40, true, false, 1e-8)
%     simulate an oscillating voltage at 0.1 Hz and 2 mV of half peak to peak voltage amplitude,
%     20 periods and 40 time points per period, taking care of reaching a stable solution,
%     and using a starting relative tolerance of 1e-8, starting from a dark
%     solution at zero applied voltage
%     calculate also the electronic current needed by impedance simulations
%   sol_i_1S_SR_is_100mHz_2mV = ISwave_EA_single_exec(sol_i_1S_SR, 2e-3, 1e-2, 20, 40, true, false, 1e-8)
%     as above, but with 1 sun illumination starting from short circuit
%     conditions
%   asymssol_i_1S_SR_is_100mHz_2mV = ISwave_EA_single_exec(asymmetricize(ssol_i_1S_SR), 2e-3, 1e-2, 20, 40, true, false, 1e-8)
%     as above, but at starting from open circuit conditions
%
% Other m-files required: pindrift, verifyStabilization, pinana
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

% shortcut
p = asymstruct_Int.p;

% suggested to use linear time mesh (1), not exponential as usual (2),
% 'cause here we're not interested in a decay, but in a oscillation going on during the whole simulated timespan
p.tmesh_type = 1;

p.tmax = periods / freq;
% it's good to simulate exactly an integer number of periods, so that the
% last point can be used as the starting point for another simulation,
% useful for reaching stabilization
p.tpoints = 1 + tpoints_per_period * periods;

if EA
    p.calcJ = 0; % in ElectroAbsorbance simulation we don't need the electrical current value
else
    % a meaningful value of calcJ should be set from pinParams
    % in Impedance Spectroscopy we need the current at the electrodes
    p.Ana = 1;
    assert(logical(p.calcJ), [mfilename ' - p.calcJ needs to be set in the input structure in order to perform an IS simulation']);
end

% accuracy
p.RelTol = RelTol;

%% define applied voltage and fitting functions
p.JV = 2; % mode for arbitrary Vapp functions

% get Efn and Efp
[Vapp_arr, ~, ~] = pinana(asymstruct_Int);

% take current voltage, as defined in pinana
Vstart = Vapp_arr(end);

% applied voltage profile
p.Vapp_params = [Vstart, deltaV, 0, 2 * pi * freq];
p.Vapp_func = @(coeff, t) coeff(1) + coeff(2) * sin(coeff(3) + coeff(4) * t);

% functions for fitting
p.J_E_func = @(coeff, t) coeff(1) + coeff(2) * sin(coeff(3) + 2 * pi * freq * t);
p.J_E_func_tilted = @(coeff, t, tilting, t_middle) p.J_E_func(coeff, t) + tilting * (t - t_middle);
% double frequency for E squared fitting, useful for Stark spectroscopy
% (ElectroAbsorption EA) simulation
p.E2_func = @(coeff, t) coeff(1) + coeff(2) * sin(-pi/2 + coeff(3) + 4 * pi * freq * t);

%% first run
disp([mfilename ' - Int: ' num2str(p.Int) '; Vdc: ' num2str(Vapp_arr(end)) ' V; Freq: ' num2str(freq) ' Hz']);
asymstruct_ISwave = pindrift(asymstruct_Int, p);

%% run until oscillating stabilization
% verify that the final time point has the same solution as the last point
% of an oscillation in the middle of the simulated timespan, so verify that "periods" was enough
% if just the first solution is requested setting reach_stability to
% false, break after the first cycle
i = 0;
if reach_stability
    warning('off', 'pindrift:verifyStabilization'); % verifyStabilization warnings are substituted by others here
    while ~verifyStabilization(asymstruct_ISwave.sol, asymstruct_ISwave.t, floor(periods / 2) / periods)
        disp([mfilename ' - solution was not stabilized, trying again'])
        asymstruct_ISwave = pindrift(asymstruct_ISwave, p);
        i = i+1;
        % in case 10 repetitions were not enough, stop
        if i > 10
            warning('pindrift:ISwave_EA_single_exec',...
                'ISwave_EA_single_exec seems that the solution did not reach complete stabilization after %s repetitions on %s periods',...
                num2str(i), num2str(periods))
            break
        end
    end
    warning('on', 'pindrift:verifyStabilization'); % re-enable warnings
end

%------------- END OF CODE --------------
