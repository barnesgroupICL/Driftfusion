function coeff = IS_EA_ana_demodulation(t, y, fun_type, freq)
%IS_EA_ANA_DEMODULATION - Calculate phase and amplitude demodulating oscillating current data from impedance spectroscopy with oscillating voltage.
% The current profile gets multiplied by the voltage profile and,
% separately, by the voltage profile with an additional 90 degrees phase.
% The integrals of the resulting profiles are related to the phase.
% This simulates the working principle of a dual-phase demodulator often
% used by lock-in amplifiers.
% https://en.wikipedia.org/wiki/Lock-in_amplifier
% https://en.wikipedia.org/wiki/Demodulation
% https://en.wikipedia.org/wiki/In-phase_and_quadrature_components
% a phase close to 90 degrees results in very small values of
% inphase_integral which can have a wrong sign and result in a wrong phase of
% -90 degrees, in that case the simulation should be repeated with smaller
% pdepe relative tolerance RelTol
%
% Syntax:  coeff = IS_EA_ana_demodulation(t, y, fun_type, freq)
%
% Inputs:
%   T - a column with the time mesh.
%   Y - a column with the value to be fitted.
%     In case of IS: current.
%   FUN_TYPE - char array, type of function as defined in fun_gen.m used
%     for generating the applied voltage time profile
%   FREQ - numeric, frequency at which the demodulation should be run
%
% Outputs:
%   COEFF - an array with the values from the demodulation: constant
%     bias, sinusoid amplitude, phase shift.
%
% Example:
%   coeff = IS_EA_ana_demodulation(ssol_i_1S_SR_is_100mHz_2mV.t', ssol_i_1S_SR_is_100mHz_2mV.Jn, 'sin', ssol_i_1S_SR_is_100mHz_2mV.par.V_fun_arg(3))
%     demodulate the components of the oscillating current
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also doIS_EA, IS_script, IS_EA_ana_fit.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

% if the input is provided in two different directions, weird output is produced
assert(size(t, 2) == 1, [mfilename ' - ' inputname(1) ' has to be provided as a single column']);
assert(size(y, 1) == size(t, 1), [mfilename ' - ' inputname(1) ' and ' inputname(2) ' need to have the same number of rows, indeed the size of the first is ' num2str(size(t)) ' while the size of the second is ' num2str(size(y))]);
% verify if we are analysing a serie of whole periods
assert(ismembertol(mod(t(end)-t(1), 1/freq),0), [mfilename ' - It is suggested to provide a whole period, including the first and last repeated points']);

func = fun_gen(fun_type);
func_fixedfreq = @(coeff,t) func([coeff(1), coeff(2), freq, coeff(3)], t);

%% profile oscillating in-phase with the voltage

% bias of applied voltage, not wanted here
fun_arg_biasV = 0;
% deltaV of applied voltage, set to one to normalize
fun_arg_deltaV = 1;
% phase in-phase with the voltage, we can use zero as it is the value set in doIS_EA
fun_arg_phase = 0; 

sin_array = func_fixedfreq([fun_arg_biasV, fun_arg_deltaV, fun_arg_phase], t);

%% profile oscillating out-of-phase with the voltage

fun_arg_phase = pi/2; % phase out-of-phase with the voltage
cos_array = func_fixedfreq([fun_arg_biasV, fun_arg_deltaV, fun_arg_phase], t);

%% modulation

% in case the input values are not oscillating around zero, subtract the
% average
% actually there should be no need to subtract the constant bias of the y input, it gets eliminated
% in the integration as constant term
max_y = max(y);
min_y = min(y);
mean_y = (min_y + max_y) / 2;
y_nobias = y - mean_y;

% for reducing the probability of numerical errors due to values close to
% the double or single var type, normalize the y
scale = 1 / (max_y - min_y);
y_rescale = y_nobias * scale;

% multiply the input profile with a purely in-phase and a purely
% out-of-phase ones
inphase = sin_array .* y_rescale;
quadrature = cos_array .* y_rescale;

%% integration

% integer the resulting profiles over an integer number of periods, if the
% integration considers an incomplete period, the result will be heavily
% affected
inphase_integral = trapz(t, inphase) / (t(end) - t(1));
quadrature_integral = trapz(t, quadrature) / (t(end) - t(1));

%% parameters extraction

% in case a confirmation of "amp" was needed, a reference for the input amplitude should be:
% disp(mean(abs(y_nobias)) * sqrt(2))
amp = 2 * sqrt(inphase_integral.^2 + quadrature_integral.^2); % Vapp_params(4) is 2pi*frequence of applied voltage
amp = amp / scale;

% get the phase via an arc-tangent
phase = atan(quadrature_integral / inphase_integral);
% as atan makes sense between -pi/2 and pi/2, it fails when the phase is
% out of this range, this fixes:
if inphase_integral < 0
    phase = phase + pi;
end
% restrict phase to nice-to-see values between -pi and pi
phase = wrapToPi(phase);

coeff = [mean_y, amp, phase];

%------------- END OF CODE --------------
