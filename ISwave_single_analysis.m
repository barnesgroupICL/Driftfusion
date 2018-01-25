function [fit_coeff, fit_i_coeff, subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t, Ji_disp] = ISwave_single_analysis(asymstruct_ISwave, minimal_mode)
%ISWAVE_SINGLE_ANALYSIS - Calculate impedance (reactance and resistance) and phase by Impedance
% Spectroscopy (ISwave) with oscillating voltage
%
% Syntax:  ISwave_single_analysis(asymstruct_ISwave, minimal_mode)
%
% Inputs:
%   ASYMSTRUCT_ISWAVE - a struct with a solution being perturbated by an
%     oscillating voltage, as generated from ISwave_single_exec
%   MINIMAL_MODE - logical, when true graphics does not get created and
%     ISwave_subtracting_analysis does not get launched, useful when
%     launched under parallelization
%
% Example:
%   ISwave_single_analysis(ISwave_single_exec(asymmetricize(ssol_i_light, 1), 1, 2e-3, 1e6, 20, 40, true, 1e-4), false)
%     do plot
%
% Other m-files required: ISwave_subtracting_analysis
% Subfunctions: none
% MAT-files required: none
%
% See also ISwave_full_exec, ISwave_single_exec, ISwave_subtracting_analysis.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

% evil shortcut
s = asymstruct_ISwave;

% round should not be needed here
% s.params.Vapp_params(4) should be pulsatance
periods = round(s.params.tmax * s.params.Vapp_params(4) / (2 * pi));
% this works if, as now is, the frequency is constant internally each solution
% round should not be needed here
fit_t_index = round((s.params.tpoints - 1) * round(periods * 0.5) / periods);
fit_t = s.t(fit_t_index:end);
fit_J = s.Jtotr(fit_t_index:end) / 1000; % in Ampere
fit_coeff = ISwave_single_fit(fit_t, fit_J, s.params.J_E_func);

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
    fit_i_coeff = ISwave_single_fit(fit_t, fit_Ji, s.params.J_E_func);
else
    fit_i_coeff = [NaN, NaN, NaN];
    Ji_disp = NaN;
end

%% plot solutions
if ~minimal_mode % disable all this stuff if under parallelization
    
    [subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t] = ISwave_subtracting_analysis(asymstruct_ISwave);

    Vapp = s.params.Vapp_func(s.params.Vapp_params, s.t);

    % fourth value of Vapp_params is pulsatance
    figure('Name', ['Single ISwave, Int ' num2str(s.params.Int) ' Freq ' num2str(s.params.Vapp_params(4) / (2 * pi))], 'NumberTitle', 'off');
        yyaxis right
        hold off
        plot(s.t, Vapp);
        ylabel('Applied voltage [V]');

        yyaxis left
        hold off
        plot(s.t, s.Jtotr / 1000); % Ampere
        hold on
        plot(s.t, s.Jdispr / 1000); % Ampere
        plot(s.t(2:end), subtracting_n_intr_t);
        plot(s.t(2:end), subtracting_n_contacts_t);
        plot(fit_t, s.params.J_E_func(fit_coeff, fit_t), 'g')
        legend_array = ["Total J", "Displacement J", "Charges intrinsic", "Charges contacts", "Fit"];
        if s.params.mui % if there was ion mobility, current due to ions have been calculated, plot stuff
            plot(s.t(2:end), subtracting_i_abs_t);
            plot(s.t(2:end), subtracting_i_t);
            plot(s.t, Ji_disp);
            plot(fit_t, s.params.J_E_func(fit_i_coeff, fit_t), 'g--');
            legend_array = [legend_array, "Displaced ions abs", "Displaced ions", "Disp J Ions", "Disp J Ions fit"];
            if max(s.params.Jpoints) % if the ionic drift has been calculated (that value is different from zero)
                plot(s.t, s.Jidrift_points(:,end) / 1000, 'k.-'); % in Ampere
                legend_array = [legend_array, "Ionic drift"];
            end
        end
        ylabel('Current [A/cm^2]');
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
