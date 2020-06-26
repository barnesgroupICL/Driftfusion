function [Voc, deltaV, monoexp_lin, monoexp_tau, biexp_lin_fast, biexp_tau_fast, biexp_lin_slow, biexp_tau_slow] = TPV_single_analysis_noMultiStart(symstruct_pulse, biexp_initial_guess_lin_fast, biexp_initial_guess_tau_fast, biexp_initial_guess_lin_slow, biexp_initial_guess_tau_slow)
%TPV_SINGLE_ANALYSIS - Fit voltage decays from single TPV solutions and do graphics
% Fit the voltage versus time profile of the provided perturbed
% solution with a mono-exponential and a bi-exponential formula.
% Try various starting point for helping the convergence of bi-exponential
% fitting.
% An initial guess from a previous fitting can be provided in input.
%
% Syntax:  [Voc, delta_v, monoexp_lin, monoexp_tau, biexp_lin_fast, biexp_tau_fast, biexp_lin_slow, biexp_tau_slow, Vslope] = TPV_single_analysis(symstruct_pulse)
%
% Inputs:
%   SYMSTRUCT_PULSE - a struct with a solution being perturbated by a light
%     pulse
%   biexp_initial_guess_lin_fast - optional, initial guess for linear
%     factor of fast decay of bi-exponential fit
%   biexp_initial_guess_tau_fast - optional, initial guess for life-time
%     of fast decay of bi-exponential fit
%   biexp_initial_guess_lin_slow - optional, initial guess for linear
%     factor of slow decay of bi-exponential fit
%   biexp_initial_guess_tau_slow - optional, initial guess for life-time
%     of slow decay of bi-exponential fit
%
% Outputs:
%   Voc - steady state voltage being present before the pulse;
%   deltaV - maximum voltage variation due to the pulse;
%   monoexp_lin - linear factor of mono-exponential fit;
%   monoexp_tau - life-time of mono-exponential fit;
%   biexp_lin_fast - linear factor of fast decay of bi-exponential fit;
%   biexp_tau_fast - life-time of fast decay of bi-exponential fit;
%   biexp_lin_slow - linear factor of slow decay of bi-exponential fit;
%   biexp_tau_slow - life-time of slow decay of bi-exponential fit.
%
% Example:
%   TPV_single_analysis(TPV_single_exec(ssol_i_light, 1e-3, 5))
%     fit and plot the output of TPV_single_exec
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also TPV_single_exec, TPVvariab_full_exec, TPVconst_full_exec.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

% increase graphics font size
set(0, 'defaultAxesFontSize', 24);
% set image dimension
set(0, 'defaultfigureposition', [0, 0, 1000, 750]);
% set line thickness
set(0, 'defaultLineLineWidth', 2);

Voc = symstruct_pulse.Voc(1);
V_minusBaseline = symstruct_pulse.Voc - Voc;

% in order to fit the decay, subset taking points after pulse off
% take value of peak of voltage caused by the pulse
% abs is included for avoiding a deltaV = 0 in cases of negative peaks
[deltaV, peakIndex] = max(abs(V_minusBaseline));
V_afterPeak = transpose(V_minusBaseline(peakIndex:end));
%Voc_during_pulse = transpose(symstruct_pulse.Voc(symstruct_pulse.t < 0)); % for slope fitting
t_afterPeak = symstruct_pulse.t(peakIndex:end);
%t_during_pulse = symstruct_pulse.t(symstruct_pulse.t < 0); % for slope fitting

% just using fit(x,y,'exp1') does not work
% anymore with matlab R2017b, now fitnlm has to be used

%% do the monoexponential fit
F = @(coef,xdata) coef(1) * exp(-xdata / coef(2)); % mono exponential function
[~, monoexp_initial_guess_tau_index] = min(abs( V_afterPeak - (deltaV / exp(1)) )); % try to guess the tau
monoexp_initial_guess = [0.9 * deltaV, t_afterPeak(monoexp_initial_guess_tau_index)]; % give approximate values for starting the fit
fit_monoexp = fitnlm(t_afterPeak', V_afterPeak, F, monoexp_initial_guess); % fitting
fit_monoexp_coef = fit_monoexp.Coefficients.Estimate; % store the fit results
monoexp_lin = fit_monoexp_coef(1);
monoexp_tau = fit_monoexp_coef(2);
monoexp_y_predicted = predict(fit_monoexp, t_afterPeak'); % generate y data for plotting the fit result

%% do the biexponential fit
% time constants are squared for avoiding negative values, as it wouldn't
% make sense otherwise
%
% amplitudes are squared for avoiding negative values as the fitting could
% easily fail. Cases where a negative peak is studied will need a code
% change
F = @(coef,xdata) coef(1)^2 * exp(-xdata / coef(2)^2) + coef(3)^2 * exp(-xdata / coef(4)^2); % bi exponential function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All this part could be rewritten using MultiStart from
% Global Optimization Toolbox
% https://mathworks.com/help/releases/R2017b/gads/multistart.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biexp_success = false;

% try biexponential fit with provided initial guess
if exist('biexp_initial_guess_lin_fast', 'var') && biexp_initial_guess_lin_fast ~= 0
    try
        biexp_initial_guess = [sqrt(biexp_initial_guess_lin_fast),...
            sqrt(biexp_initial_guess_tau_fast), sqrt(biexp_initial_guess_lin_slow),...
            sqrt(biexp_initial_guess_tau_slow)]; % give approximate values for starting the fit
        fit_biexp = fitnlm(t_afterPeak', V_afterPeak, F, biexp_initial_guess); % fitting
        biexp_success = true;
    catch % in case biexp fit fails no problem
        warning('Biexponential fit failed for light intensity %s', num2str(symstruct_pulse.p.Int));
    end
end

% if previous fit failed, try biexponential fit with meaningful initial guess
if ~biexp_success
    %try
        % try to guess the slow component parameters
        % slow component intensity is usually between 10% and 30% of the toal peak
        [~, biexp_initial_guess_tau_slow_index] = min(abs( V_afterPeak - (deltaV * 0.2 * exp(-5)) ));
        biexp_initial_guess_tau_slow = t_afterPeak(biexp_initial_guess_tau_slow_index) * 10;
        biexp_initial_guess_lin_slow = V_afterPeak(biexp_initial_guess_tau_slow_index) * exp(5);
        % try to guess the fast component parameters
        V_afterPeak_minusSlow = V_afterPeak - biexp_initial_guess_lin_slow * exp(-t_afterPeak / biexp_initial_guess_tau_slow);
        deltaV_minusSlow = max(V_afterPeak_minusSlow);
        [~, biexp_initial_guess_tau_fast_index] = min(abs(V_afterPeak_minusSlow - (deltaV_minusSlow * exp(-1/5)) ));
        biexp_initial_guess_tau_fast = t_afterPeak(biexp_initial_guess_tau_fast_index);
        biexp_initial_guess_lin_fast = V_afterPeak_minusSlow(biexp_initial_guess_tau_fast_index) * exp(1/5);
        biexp_initial_guess = [sqrt(biexp_initial_guess_lin_fast),...
            sqrt(biexp_initial_guess_tau_fast), sqrt(biexp_initial_guess_lin_slow),...
            sqrt(biexp_initial_guess_tau_slow)]; % give approximate values for starting the fit
    % for debugging the initial guess
    %    format short e
    %    disp(biexp_initial_guess);
    %    figure(200)
    %    hold off
    %    plot(symstruct_pulse.t, V_minusBaseline, '.')
    %    hold on
    %    fplot(@(t) biexp_initial_guess_lin_fast * exp( -t / biexp_initial_guess_tau_fast ) + biexp_initial_guess_lin_slow * exp( -t / biexp_initial_guess_tau_slow ), [0, symstruct_pulse.t(end)]);
    %    hold off
    %    format short
        fit_biexp = fitnlm(t_afterPeak', V_afterPeak, F, biexp_initial_guess); % fitting
        biexp_success = true;
    %catch % in case biexp fit fails no problem
    %    warning('Biexponential fit failed for light intensity %s', num2str(symstruct_pulse.p.Int));
    %end
end

% if previous fit failed, try biexponential fit with initial guess based on monoexponential fit
if ~biexp_success
    try
        % give some values for starting the fit
        biexp_initial_guess = [sqrt(abs(monoexp_lin) / 2), sqrt(abs(monoexp_tau) * 10), sqrt(abs(monoexp_lin) / 2), sqrt(abs(monoexp_tau) / 10)];
        fit_biexp = fitnlm(t_afterPeak', V_afterPeak, F, biexp_initial_guess); % fitting
        biexp_success = true;
    catch % in case biexp fit fails again
        warning('Biexponential fit failed for light intensity %s', num2str(symstruct_pulse.p.Int));
    end
end

if biexp_success
    fit_biexp_coef = fit_biexp.Coefficients.Estimate; % store the fit results
    % get fast decay
    [biexp_tau_fast, fast_coef_index] = min(fit_biexp_coef([2,4]).^2);
    biexp_lin_fast = fit_biexp_coef(fast_coef_index * 2 - 1)^2;
    [biexp_tau_slow, slow_coef_index] = max(fit_biexp_coef([2,4]).^2);
    biexp_lin_slow = fit_biexp_coef(slow_coef_index * 2 - 1)^2;
    biexp_y_predicted = predict(fit_biexp, t_afterPeak'); % generate y data for plotting the fit result
else
    % just plot a point at NaN (in case of graphs it just disappears)
    biexp_lin_fast = 0;
    biexp_tau_fast = 0;
    biexp_lin_slow = 0;
    biexp_tau_slow = 0;
    biexp_y_predicted = 0;
end

%% calculate the slope during the pulse
% 
% % fit with linear
% fit_slope = fitlm(t_during_pulse', Voc_during_pulse);
% fit_slope_coef = fit_slope.Coefficients.Estimate;
% %disp(fit_slope_coef)
% slope_initial_guess = [fit_slope_coef(2) / 1e-5, 1e-5];
% start_time = symstruct_pulse.t(1);
% 
% % Fslope = @(coef,xdata)coef(1)*(1-exp(-(xdata+start_time)/coef(2))); % exponential function
% % try
% %     fit_slope = fitnlm(t_during_pulse', Voc_during_pulse, Fslope, slope_initial_guess);
% % catch % in case exp fit fails
% %     warning('Slope exponential fit failed for light intensity %s', num2str(symstruct_pulse.p.Int));
% % end
% timestep = symstruct_pulse.t(2) - symstruct_pulse.t(1);
% 
% Vstep = predict(fit_slope, symstruct_pulse.t(2)) - predict(fit_slope, symstruct_pulse.t(1));
% Vslope = Vstep / timestep;
% 
% Vslope_y_predicted = predict(fit_slope, t_during_pulse'); % generate y data for plotting the fit result


%% plot figures

figure('Name', ['Single TPV at light intensity ' num2str(symstruct_pulse.p.Int)], 'NumberTitle', 'off');
    hold off
    semilogy(t_afterPeak, monoexp_y_predicted, 'r');
    hold on
    semilogy(symstruct_pulse.t, V_minusBaseline, '.');
    semilogy(t_afterPeak, biexp_y_predicted, 'g');
    [~, lim_point_index] = min(abs(t_afterPeak - 50 * biexp_tau_fast)); % fifty times more the Tau of fast decay
    % verifies that biexp fit was succesful && if the y limit is going to
    % be valid && if the resulting fit approaches zero at the end of times
    % otherwise use automatic xlim and ylim
    if biexp_tau_fast && V_afterPeak(lim_point_index) > 0 && F(fit_biexp_coef, t_afterPeak(end)) < deltaV / 20
        % set time limits, shorter at the end for making more clear the fast decay
        xlim([symstruct_pulse.t(1), t_afterPeak(lim_point_index)]);
        ylim([V_afterPeak(lim_point_index) / 1.5, deltaV * 1.5]); % set Voc limits, ignore tail
    else
        xlim([symstruct_pulse.t(1), symstruct_pulse.t(end)]);
        ylim([min(V_afterPeak(V_afterPeak > 0)) / 1.5, deltaV * 1.5]);
    end
    
    % debug linear fit
    %fplot(@(t) fit_slope_coef(1) + fit_slope_coef(2) * t, [start_time, 0]);
    
    % debug initial guess for slope fit
    %fplot(@(t) slope_initial_guess(1)*(1-exp(-(t + start_time)/slope_initial_guess(2) )), [start_time, 0]);

    % plot the slope
    %semilogy(t_during_pulse, Vslope_y_predicted, 'g');
    
    xlabel('Time [s]');
    ylabel('Voltage relative to Voc [V]');
    legend('mono exp fit', 'TPV', 'bi exp fit');
    hold off
    figure('Name', 'Stocazzo', 'NumberTitle', 'off')
    
%------------- END OF CODE --------------
