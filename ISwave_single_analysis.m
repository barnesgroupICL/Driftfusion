function [coeff, i_coeff, subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t, Ji_disp] = ISwave_single_analysis(asymstruct_ISwave, minimal_mode, demodulation)
%ISWAVE_SINGLE_ANALYSIS - Calculate impedance (reactance and resistance) and phase by Impedance
% Spectroscopy (ISwave) with oscillating voltage
%
% Syntax:  ISwave_single_analysis(asymstruct_ISwave, minimal_mode, demodulation)
%
% Inputs:
%   ASYMSTRUCT_ISWAVE - a struct with a solution being perturbated by an
%     oscillating voltage, as generated from ISwave_EA_single_exec
%   MINIMAL_MODE - logical, when true graphics does not get created and
%     ISwave_subtracting_analysis does not get launched, useful when
%     launched under parallelization
%   DEMODULATION - logical, get phase via demodulation instead of using a fitting
%
% Outputs:
%  
% Example:
%   ISwave_single_analysis(ISwave_EA_single_exec(asymmetricize(ssol_i_light, 1), 1, 2e-3, 1e6, 20, 40, true, false, 1e-4), false, true)
%     do plot
%
% Other m-files required: ISwave_subtracting_analysis, ISwave_EA_single_fit, ISwave_EA_single_demodulation
% Subfunctions: none
% MAT-files required: none
%
% See also ISwave_full_exec, ISwave_EA_single_exec, ISwave_subtracting_analysis, ISwave_EA_single_demodulation, ISwave_EA_single_fit.

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

% evil shortcut
s = asymstruct_ISwave;

% round should not be needed here
% s.p.Vapp_params(4) should be pulsatance
periods = round(s.p.tmax * s.p.Vapp_params(4) / (2 * pi));

% here is critical that exactely an entire number of periods is provided to
% ISwave_single_demodulator
fit_t_index = round((s.p.tpoints - 1) * round(periods / 2) / periods) + 1;
fit_t = s.t(fit_t_index:end)';
fit_J = s.Jn(fit_t_index:end) / 1000; % in Ampere

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
    coeff = ISwave_EA_single_demodulation(fit_t, fit_J_flat, s.p.Vapp_func, s.p.Vapp_params);
else
    coeff = ISwave_EA_single_fit(fit_t, fit_J_flat, s.p.J_E_func);
end

%% calculate ionic contribution
if s.p.mui % if there was ion mobility, current due to ions have been calculated, fit it
    % Intrinsic points logical array
    itype_points= (s.x >= s.p.tp & s.x <= s.p.tp + s.p.ti);
    % subtract background ion concentration, for having less noise in trapz
    i_matrix = s.sol(:, :, 3) - s.p.NI;
    % calculate electric field due to ions
    Efield_i = s.p.e * cumtrapz(s.x, i_matrix, 2) / s.p.eppi;
    % an average would be enough if the spatial mesh was homogeneous in the
    % intrinsic, indeed I have to use trapz for considering the spatial mesh
    Efield_i_mean = trapz(s.x(itype_points), Efield_i(:, itype_points), 2) / s.p.ti;
    % calculate displacement current due to ions
    Ji_disp = -s.p.eppi * gradient(Efield_i_mean, s.t); % in Amperes

    fit_Ji = Ji_disp(fit_t_index:end); % in Ampere
    
    % remove some tilting from fit_Ji to get better fit and better demodulation in case of
    % unstabilized solutions
    % this assumes that the first and last point are at the same point in the oscillating voltage
    delta_Ji = fit_Ji(end) - fit_Ji(1);
    tilting_i = delta_Ji/delta_t;
    % because of this, the bias value that will be obtained from the fit/demodulation is not going to be correct
    % here the fit_Ji has to be provided as a row
    fit_Ji_flat = fit_Ji - tilting_i * (fit_t - t_middle);

    if demodulation
        i_coeff = ISwave_EA_single_demodulation(fit_t, fit_Ji_flat, s.p.Vapp_func, s.p.Vapp_params);
    else
        i_coeff = ISwave_EA_single_fit(fit_t, fit_Ji_flat, s.p.J_E_func);
    end
else
    i_coeff = [NaN, NaN, NaN];
    Ji_disp = NaN;
end

%% plot solutions
if ~minimal_mode % disable all this stuff if under parallelization or if explicitly asked to not plot any graphics
    
    [subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t,...
        subtracting_i_abs_t, subtracting_i_t] = ISwave_subtracting_analysis(asymstruct_ISwave);

    Vapp = s.p.Vapp_func(s.p.Vapp_params, s.t);

    % fourth value of Vapp_params is pulsatance
    figure('Name', ['Single ISwave, Int ' num2str(s.p.Int) ' Freq '...
        num2str(s.p.Vapp_params(4) / (2 * pi))], 'NumberTitle', 'off');
        yyaxis right
        hold off
        i=1;
        h(i) = plot(s.t, Vapp, 'r', 'LineWidth', 2);
        legend_array = "Applied Voltage";
        ylabel('Applied voltage [V]');

        yyaxis left
        hold off
        i=i+1; h(i) = plot(s.t, s.Jn, 'b--', 'LineWidth', 2); % mA
        hold on
        %i=i+1; h(i) = plot(s.t, -s.Jdispr); % mA
        i=i+1; h(i) = plot(s.t(2:end), -subtracting_n_intr_t * 1000); % mA
        i=i+1; h(i) = plot(s.t(2:end), -subtracting_n_contacts_t * 1000); % mA
        i=i+1; h(i) = plot(fit_t, s.p.J_E_func_tilted(coeff, fit_t, tilting, t_middle) * 1000, 'kx-'); % mA

        legend_array = [legend_array, "Current", "Charge variation intrinsic", "Charge variation contacts", "1st harmonic"]; % "Displacement current"
        if s.p.mui % if there was ion mobility, current due to ions have been calculated, plot stuff
            i=i+1; h(i) = plot(s.t, Ji_disp * 1000); % mA
            i=i+1; h(i) = plot(fit_t, s.p.J_E_func_tilted(i_coeff, fit_t, tilting_i, t_middle) * 1000, 'g--'); % mA
            legend_array = [legend_array, "Ionic displacement current", "Ionic disp J fit"];
        end
        ylabel('Current [mA/cm^2]');
        hold off
        xlabel('Time [s]');
        legend(h, legend_array);
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
