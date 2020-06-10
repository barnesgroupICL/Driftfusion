function coeff = IS_EA_ana_fit(t, y, fun_type, freq)
%IS_EA_ANA_FIT - Calculate phase and amplitude fitting oscillating current data from impedance spectroscopy with oscillating voltage
% This is an alternative to IS_EA_ana_demodulation.m which can be used to confirm the demodulation results.
%
% Syntax:  fit_coeff = IS_EA_ana_fit(t, y, fun_type, freq)
%
% Inputs:
%   T - an array with the time mesh.
%   Y - an array with the value to be fitted.
%     In case of IS: current.
%   FUN_TYPE - an anonymous function to be used for the fitting, containing
%     a constant bias, the amplitude of the sinusoid and its phase. The
%     frequency is taken as it is from input, it is not fitted.
%   FREQ - numeric, frequency at which the fit function will be
%
% Outputs:
%   COEFF - an array with the values from the sinusoidal fit: constant
%     bias, sinusoid amplitude, phase shift.
%
% Example:
%   coeff = IS_EA_ana_fit(ssol_i_1S_SR_is_100mHz_2mV.t, ssol_i_1S_SR_is_100mHz_2mV.Jn, 'sin', ssol_i_1S_SR_is_100mHz_2mV.par.V_fun_arg(3))
%     fits the charge vs voltage data
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also doIS_EA, IS_script, IS_EA_ana_demodulation.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: December 2017

%------------- BEGIN CODE --------------

% check input
assert(size(t, 1) == 1 || size(t, 2) == 1, [mfilename ' - ' inputname(1) ' has to be provided as a single column or row'])
assert(size(y, 1) == 1 || size(y, 2) == 1, [mfilename ' - ' inputname(2) ' has to be provided as a single column or row'])

% solution taken from https://mathworks.com/matlabcentral/answers/150847-problem-using-nlinfit-poor-fit#answer_148643
% basically the fit can be bad because we're close to the single floating point
% precision limit, rescaling the values is enough for having a good fit
scale = 1 / max(abs(y));
y_rescale = y * scale;

func = fun_gen(fun_type);
func_fixedfreq = @(coeff,t) func([coeff(1), coeff(2), freq, coeff(3)], t);

%% initial guess
ig_bias = mean(y_rescale);
ig_amp = range(y_rescale) / 2;
ig_phase = pi / 2;

initial_guess = [ig_bias, ig_amp, ig_phase];

%% fitting complete model
% non linear fit
fit = fitnlm(t, y_rescale, func_fixedfreq, initial_guess);
bias = fit.Coefficients.Estimate(1) / scale;
amp = fit.Coefficients.Estimate(2) / scale;
phase = fit.Coefficients.Estimate(3);

% a negative amplitude is ugly, invert it and shift phase
if amp < 0
    amp = -1 * amp;
    phase = phase + pi;
end

% obtain nice values of phase
phase = wrapToPi(phase);

coeff = [bias, amp, phase];

%------------- END OF CODE --------------
