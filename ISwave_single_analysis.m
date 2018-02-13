function [coeff, i_coeff, subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t, Ji_disp] = ISwave_single_analysis(asymstruct_ISwave, minimal_mode, demodulation)
%ISWAVE_SINGLE_ANALYSIS - Calculate impedance (reactance and resistance) and phase by Impedance
% Spectroscopy (ISwave) with oscillating voltage
%
% Syntax:  ISwave_single_analysis(asymstruct_ISwave, minimal_mode, demodulation)
%
% Inputs:
%   ASYMSTRUCT_ISWAVE - a struct with a solution being perturbated by an
%     oscillating voltage, as generated from ISwave_single_exec
%   MINIMAL_MODE - logical, when true graphics does not get created and
%     ISwave_subtracting_analysis does not get launched, useful when
%     launched under parallelization
%   DEMODULATION - logical, get phase via demodulation instead of using a fitting
%
% Example:
%   ISwave_single_analysis(ISwave_single_exec(asymmetricize(ssol_i_light, 1), 1, 2e-3, 1e6, 20, 40, true, true, 1e-4), false, true)
%     do plot
%
% Other m-files required: ISwave_subtracting_analysis, ISwave_single_fit, ISwave_single_demodulation
% Subfunctions: none
% MAT-files required: none
%
% See also ISwave_full_exec, ISwave_single_exec, ISwave_subtracting_analysis, ISwave_single_demodulation, ISwave_single_fit.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

% increase graphics font size
set(0, 'defaultAxesFontSize', 30);
% set image dimension
set(0, 'defaultfigureposition', [0, 0, 1000, 750]);
% set line thickness
set(0, 'defaultLineLineWidth', 2);

% evil shortcut
s = asymstruct_ISwave;

% round should not be needed here
% s.params.Vapp_params(4) should be pulsatance
periods = round(s.params.tmax * s.params.Vapp_params(4) / (2 * pi));
% this works if, as now is, the frequency is constant internally each solution
% round should not be needed here
% here is critical that exactely an entire number of periods is provided to
% ISwave_single_demodulator
fit_t_index = round((s.params.tpoints - 1) * round(periods * 0.5) / periods) + 1;
fit_t = s.t(fit_t_index:end);
fit_J = s.Jtotr(fit_t_index:end) / 1000; % in Ampere

% remove some tilting from fit_J to get better fit and better demodulation in case of
% unstabilized solutions. In case of noisy solutions this could work badly
% but should not affect too much the fitting/demodulation
delta_t = fit_t(end) - fit_t(1);
% this assumes that the first and last point are at the same point in the oscillating voltage
delta_J = fit_J(end) - fit_J(1);
tilting = delta_J/delta_t;
t_middle = fit_t(round(end/2));
% because of this, the bias value that will be obtained from the fit/demodulation is not going to be correct
fit_J_flat = fit_J - tilting * (fit_t - t_middle);

if demodulation
    coeff = ISwave_single_demodulation(fit_t, fit_J_flat, s.params.Vapp_func, s.params.Vapp_params);
else
    coeff = ISwave_single_fit(fit_t, fit_J_flat, s.params.J_E_func);
end

%% calculate ionic contribution
if s.params.mui % if there was ion mobility, current due to ions have been calculated, fit it
    % Intrinsic points logical array
    itype_points= (s.x >= s.params.tp & s.x <= s.params.tp + s.params.ti);
    % subtract background ion concentration, for having less noise in trapz
    i_matrix = s.sol(:, :, 3) - s.params.NI;
    % calculate electric field due to ions
    Efield_i = s.params.e * cumtrapz(s.x, i_matrix, 2) / s.params.eppi;
    % an average would be enough if the spatial mesh was homogeneous in the
    % intrinsic, indeed I have to use trapz for considering the spatial mesh
    Efield_i_mean = trapz(s.x(itype_points), Efield_i(:, itype_points), 2) / s.params.ti;
    % calculate displacement current due to ions
    Ji_disp = s.params.eppi * gradient(Efield_i_mean, s.t); % in Amperes

    fit_Ji = Ji_disp(fit_t_index:end); % in Ampere
    
    % remove some tilting from fit_Ji to get better fit and better demodulation in case of
    % unstabilized solutions
    % this assumes that the first and last point are at the same point in the oscillating voltage
    delta_Ji = fit_Ji(end) - fit_Ji(1);
    tilting_i = delta_Ji/delta_t;
    % because of this, the bias value that will be obtained from the fit/demodulation is not going to be correct
    % here the fit_Ji has to be provided as a row
    fit_Ji_flat = fit_Ji' - tilting_i * (fit_t - t_middle);

    if demodulation
        i_coeff = ISwave_single_demodulation(fit_t, fit_Ji_flat, s.params.Vapp_func, s.params.Vapp_params);
    else
        i_coeff = ISwave_single_fit(fit_t, fit_Ji_flat, s.params.J_E_func);
    end
else
    i_coeff = [NaN, NaN, NaN];
    Ji_disp = NaN;
end

%% plot solutions
if ~minimal_mode % disable all this stuff if under parallelization or if explicitly asked to not plot any graphics
    
    [subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t,...
        subtracting_i_abs_t, subtracting_i_t] = ISwave_subtracting_analysis(asymstruct_ISwave);

    Vapp = s.params.Vapp_func(s.params.Vapp_params, s.t);

    % fourth value of Vapp_params is pulsatance
    figure('Name', ['Single ISwave, Int ' num2str(s.params.Int) ' Freq '...
        num2str(s.params.Vapp_params(4) / (2 * pi))], 'NumberTitle', 'off');
        yyaxis right
        hold off
        plot(s.t, Vapp, 'r', 'LineWidth', 2);
        ylabel('Applied voltage [V]');

        yyaxis left
        hold off
        plot(s.t, -s.Jtotr, 'LineWidth', 2); % mA
        hold on
        plot(s.t, -s.Jdispr); % mA
        plot(s.t(2:end), -subtracting_n_intr_t * 1000); % mA
        plot(s.t(2:end), -subtracting_n_contacts_t * 1000); % mA
        plot(fit_t, -s.params.J_E_func_tilted(coeff, fit_t, tilting, t_middle) * 1000, 'k') % mA

        legend_array = ["Current", "Displacement J", "Charges intrinsic", "Charges contacts", "1st harmonic"];
        if s.params.mui % if there was ion mobility, current due to ions have been calculated, plot stuff
            % plot(s.t(2:end), -subtracting_i_abs_t);
            % plot(s.t(2:end), -subtracting_i_t);
            plot(s.t, -Ji_disp * 1000); % mA
            plot(fit_t, -s.params.J_E_func_tilted(i_coeff, fit_t, tilting_i, t_middle) * 1000, 'g--'); % mA
            legend_array = [legend_array, "Disp J Ions", "J ions fit"];
            % if max(s.params.Jpoints) % if the ionic drift has been calculated (that value is different from zero)
                % plot(s.t, -s.Jidrift_points(:,end) / 1000, 'k.-'); % in Ampere
                % legend_array = [legend_array, "Ionic drift"];
            % end
        end
        ylabel('Current [mA/cm^2]');
        hold off
        xlabel('Time [s]');
        legend(legend_array);
        hold off
else
    subtracting_n_t = NaN;
    subtracting_n_intr_t = NaN;
    subtracting_n_contacts_t = NaN;
    subtracting_i_abs_t = NaN;
    subtracting_i_t = NaN;
    Ji_disp = NaN;
end

%------------- END OF CODE --------------
