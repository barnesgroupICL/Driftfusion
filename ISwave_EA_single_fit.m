function fit_coeff = ISwave_single_fit(t, y, J_E_func)
%ISWAVE_SINGLE_FIT - Calculate phase and amplitude fitting oscillating current data from impedance spectroscopy with oscillating voltage
%
% Syntax:  fit_coeff = ISwave_single_fit(t, y, func)
%
% Inputs:
%   T - an array with the time mesh.
%   Y - an array with the value to be fitted.
%     In case of IS: current in Ampere.
%   J_E_FUNC - fitting function containing a constant bias, the amplitude of
%     the sinusoid and its phase. The frequency is taken as it is from
%     input, it is not fitted.
%
% Outputs:
%   FIT_COEFF - an array with the values from the sinusoidal fit: constant
%     bias, sinusoid amplitude, phase shift.
%
% Example:
%   ISwave_single_fit(ssol_i_light_Int_1_Freq_100_ISwave.t, ssol_i_light_Int_1_Freq_100_ISwave.Jtotr, @(coeff, t) coeff(1) + coeff(2) * sin(coeff(3) + 2 * pi * 100 * t))
%     fits the charge vs voltage data
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also ISwave_single_exec, ISwave_full_exec, ISwave_single_demodulation.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: December 2017

%------------- BEGIN CODE --------------

assert(size(t, 1) == 1 || size(t, 2) == 1, [mfilename ' - ' inputname(1) ' has to be provided as a single column or row'])
assert(size(y, 1) == 1 || size(y, 2) == 1, [mfilename ' - ' inputname(2) ' has to be provided as a single column or row'])

% solution taken from https://mathworks.com/matlabcentral/answers/150847-problem-using-nlinfit-poor-fit#answer_148643
% basically the fit was bad 'cause we're close to the single floating point
% precision limit, rescaling the values is enough for fixing
scale = 1 / max(abs(y));
y_rescale = y * scale;

%% starting conditions
ig_bias = mean(y_rescale);
ig_amp = range(y_rescale) / 2;
ig_phase = pi / 2;

initial_guess = [ig_bias, ig_amp, ig_phase];

% for printing graphics about initial guesses
%     disp(['scale ' num2str(scale)])
%     figure('Name', 'Initial Guess', 'NumberTitle', 'off')
%     hold off
%     plot(t, y, '+')
%     hold on
%     plot(t, J_E_func(initial_guess, t) / scale)
%    legend(["data", "initial guess"])

%% fitting complete model
fit = fitnlm(t, y_rescale, J_E_func, initial_guess);
bias = fit.Coefficients.Estimate(1) / scale;
amp = fit.Coefficients.Estimate(2) / scale;
phase = fit.Coefficients.Estimate(3);
if amp < 0
    amp = -1 * amp;
    phase = phase + pi;
end
phase = wrapToPi(phase);

fit_coeff = [bias, amp, phase];

% fitting complete model
%     figure('Name', 'Fit Results', 'NumberTitle', 'off')
%     hold off
%     plot(t, y, '+')
%     hold on
%     plot(t, predict(fit, transpose(t)) / scale)
%    legend(["data", "fit"])


%------------- END OF CODE --------------
