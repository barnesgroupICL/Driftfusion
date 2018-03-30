function coeff = ISwave_EA_single_demodulation(t, y, Vapp_func, Vapp_params)
%ISWAVE_SINGLE_DEMODULATION - Calculate phase and amplitude demodulating
% oscillating current data from impedance spectroscopy with oscillating
% voltage. The current profile gets multiplied by the voltage profile and,
% separately, by the voltage profile with an additional 90 degrees phase,
% the integrals of the resulting profiles are related to the phase.
% This simulates the working principle of a dual-phase demodulator often
% used by lock-in amplifiers.
% https://en.wikipedia.org/wiki/Lock-in_amplifier
% https://en.wikipedia.org/wiki/Demodulation
% https://en.wikipedia.org/wiki/In-phase_and_quadrature_components
%
% Syntax:  coeff = ISwave_EA_single_demodulation(t, y, Vapp_func, Vapp_params)
%
% Inputs:
%   T - an array with the time mesh.
%   Y - a row array with the value to be fitted.
%     In case of IS: current in Ampere.
%   VAPP_FUNC - function containing a constant bias, the amplitude of
%     the sinusoid, its phase and the frequency
%   VAPP_PARAMS - are the parameters used for calculating the applied
%     oscillating voltage with VAPP_FUNC
%
% Outputs:
%   COEFF - an array with the values from the demodulation: constant
%     bias, sinusoid amplitude, phase shift.
%
% Example:
%   ISwave_EA_single_demodulation(ssol_i_light_Int_1_Freq_100_ISwave.t, ssol_i_light_Int_1_Freq_100_ISwave.Jtotr, @(coeff, t) coeff(1) + coeff(2) * sin(coeff(3) + coeff(4) * t), ssol_i_light_Int_1_Freq_100_ISwave.)
%     demodulate the components of the oscillating current
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also ISwave_EA_single_exec, ISwave_full_exec, ISwave_EA_single_fit.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: December 2017

%------------- BEGIN CODE --------------

% if the input is provided in two different directions, weird output is produced
assert(size(t, 2) == 1, [mfilename ' - ' inputname(1) ' has to be provided as a single column']);
assert(size(y, 1) == size(t, 1), [mfilename ' - ' inputname(1) ' and ' inputname(2) ' need to have the same number of rows']);
% verify if the first and last point are identical, or into the double
% precision error
if Vapp_func(Vapp_params, t(1)) >= Vapp_func(Vapp_params, t(end)) + 1e-14 || Vapp_func(Vapp_params, t(1)) <= Vapp_func(Vapp_params, t(end)) - 1e-14
    warning('pindrift:demodulation', 'It is suggested to provide a whole period, including the first and last repeated points')
end

Vapp_params(1) = 0; % Vapp_params(1) is bias of applied voltage, not wanted here
Vapp_params(2) = 1; % Vapp_params(2) is deltaV of applied voltage, set to one to normalize

sin_array = Vapp_func(Vapp_params, t);
cos_params = Vapp_params;
cos_params(3) = cos_params(3) + pi/2; % Vapp_params(3) is phase of applied voltage
cos_array = Vapp_func(cos_params, t);

% in case the provided values are not oscillating around zero, subtract the
% average
max_y = max(y);
min_y = min(y);
mean_y = (min_y + max_y) / 2;
y_nobias = y - mean_y;

% for reducing the probability of numerical errors due to values close to
% the double or single var type, normalize the y
scale = 1 / (max_y - min_y);
y_rescale = y_nobias * scale;

% no need to subtract the constant bias of the y input, it gets eliminated
% in the integration as constant term
inphase = sin_array .* y_rescale;
quadrature = cos_array .* y_rescale;

inphase_integral = trapz(t, inphase) / (t(end) - t(1));
quadrature_integral = trapz(t, quadrature) / (t(end) - t(1));

% a reference for the input amplitude should be:
% disp(mean(abs(y_nobias)) * sqrt(2))
amp = 2 * sqrt(inphase_integral.^2 + quadrature_integral.^2); % Vapp_params(4) is 2pi*frequence of applied voltage
amp = amp / scale;

phase = atan(quadrature_integral / inphase_integral);
% as atan makes sense between -pi/2 and pi/2, it fails when the phase is
% out of this range, this fixes:
if inphase_integral < 0
    phase = phase + pi;
end
% restrict phase to nice-to-see values between -pi and pi
phase = wrapToPi(phase);

% a phase close to 90 degrees results in very small values of
% inphase_integral
% which can have a wrong sign and result in a wrong phase of
% -90 degrees
coeff = [mean_y, amp, phase];

%------------- END OF CODE --------------
