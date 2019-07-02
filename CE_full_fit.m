function fit_coeff = CE_full_fit(delta_voltages, charges)
%CE_FULL_FIT - Calculate geometric and chemical capacitance from CE data
%
% Syntax:  fit_coeff = CE_full_fit(delta_voltages, charges)
%
% Inputs:
%   DELTA_VOLTAGES - an array containing the voltage drops which were
%     imposed in the CE experiments. As usually the CE experiments goes from
%     open circuit to short circuit, the voltage drop is equal to the VOC.
%   CHARGES - an array with the extracted charges from CE at various light
%     intensities.
%
% Outputs:
%   FIT_COEFF - an array with the values from the linear + exponential fit:
%     linear gradient, exponential function linear factor, exponent.
%
% Example:
%   CE_full_fit(CE_struct.Voc, CE_struct.extracting_charges)
%     fits the charge vs voltage data
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also CE_full_exec, CE_full_analysis.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: December 2017

%------------- BEGIN CODE --------------

% easier to have sorted values
[delta_voltages, sort_index] = sort(delta_voltages);
charges = charges(sort_index);

% solution taken from https://mathworks.com/matlabcentral/answers/150847-problem-using-nlinfit-poor-fit#answer_148643
% basically the fit was bad 'cause we're close to the single floating point
% precision limit, rescaling the values is enough for fixing
scale = 1 / charges(ceil(end/2)); % in case of just one point, takes the only value
charges_rescale = charges * scale;

disp([mfilename ' - full fit'])
%% starting conditions for exp + lin fit
lin = charges_rescale(ceil(end/2)) / delta_voltages(ceil(end/2)); % in case of just one point, takes the only value
ig_exp_gamma = 26;
ig_exp_amp = (charges_rescale(end) - lin * delta_voltages(end)) / (exp(delta_voltages(end) * ig_exp_gamma) - 1);
initial_guess = [lin, ig_exp_amp, ig_exp_gamma];

% real fitting function
F = @(coef, xdata) coef(1) * xdata + coef(2) * (exp(xdata * coef(3)) - 1);

% for printing graphics about initial guesses
%     figure(1)
%     hold off
%     plot(delta_voltages, charges_rescale, '+')
%     hold on
%     plot(delta_voltages, F(initial_guess, delta_voltages))
%     % for printing the initial guess
%     disp([num2str(initial_guess(1) / scale) ', ' num2str(initial_guess(2) / scale) ', ' num2str(initial_guess(3))]);

%% fitting complete model
fit = fitnlm(delta_voltages, charges_rescale, F, initial_guess);

% create output
fit_coeff1 = fit.Coefficients.Estimate(1) / scale;
fit_coeff2 = fit.Coefficients.Estimate(2) / scale;
fit_coeff3 = fit.Coefficients.Estimate(3);
fit_coeff = [fit_coeff1, fit_coeff2, fit_coeff3];

% for printing graphics about initial guesses
%     figure(2)
%     hold off
%     plot(delta_voltages, charges_rescale, '+')
%     hold on
%     plot(delta_voltages, predict(fit, delta_voltages))
%     % for printing the fit results
%     disp([num2str(fit_coeff1) ', ' num2str(fit_coeff2) ', ' num2str(fit_coeff3)]);

%------------- END OF CODE --------------
