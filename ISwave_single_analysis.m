function [fit_coeff, fit_i_coeff, subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t, Ji_disp] = ISwave_single_analysis(sol_i_Int_ISwave, minimal_mode)
% calculate impedance (reactance and resistance) and phase by impedance
% spectroscopy with oscillating voltage

% when run in parallelization, minimal_mode is set so that graphics does
% get created (they does not appear when under parallelization) and
% subtracting_analysis does not get launched

% evil shortcut
s = sol_i_Int_ISwave;

% round should not be needed here
% s.params.Vapp_params(4) should be pulsatance
periods = round(s.params.tmax * s.params.Vapp_params(4) / (2 * pi));
% this works if, as now is, the frequency is constant internally each solution
% round should not be needed here
fit_t_index = round((s.params.tpoints - 1) * round(periods * 0.5) / periods);
fit_t = s.t(fit_t_index:end);
fit_J = s.Jtotr(fit_t_index:end) / 1000; % in Ampere
fit_coeff = ISwave_single_fit(fit_t, fit_J, s.params.J_func);

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
    fit_i_coeff = ISwave_single_fit(fit_t, fit_Ji, s.params.J_func);
else
    fit_i_coeff = [NaN, NaN, NaN];
end

if ~minimal_mode % disable all this stuff if under parallelization
    
    [subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t] = ISwave_subtracting_analysis(sol_i_Int_ISwave);

    Vapp = s.params.Vapp_func(s.params.Vapp_params, s.t);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate current out of phase

    % calculate amplitude of out of phase oscillation
    J_amp_outphase = fit_coeff(2) * sin(fit_coeff(3));

    % calculate bias of out of phase component as fraction of total J_bias
    J_bias = fit_coeff(1);
    J_amp = fit_coeff(2);
    J_bias_outphase = J_bias * J_amp_outphase / J_amp;
    J_outphase = s.params.J_func([J_bias_outphase, J_amp_outphase, pi/2], fit_t);

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
        plot(fit_t, s.params.J_func(fit_coeff, fit_t), 'g')
        plot(fit_t, J_outphase, 'r')
        legend_array = ["Total J", "Displacement J", "Charges intrinsic", "Charges contacts", "Fit", "Out of Phase Current"];
        if s.params.mui % if there was ion mobility, current due to ions have been calculated, plot stuff
            plot(s.t(2:end), subtracting_i_abs_t);
            plot(s.t(2:end), subtracting_i_t);
            plot(s.t, Ji_disp);
            legend_array = [legend_array, "Displaced ions abs", "Displaced ions", "Disp Curr Ions"];
            if max(s.params.Jpoints) % if the ionic drift has been calculated (that value is different from zero)
                plot(s.t, s.Jidrift_points(:,end) / 1000, 'k.-'); % in Ampere
                plot(fit_t, s.params.J_func(fit_i_coeff, fit_t), 'g--');
                legend_array = [legend_array, "Ionic drift", "Ion drift fit"];
            end
        end
        ylabel('Current [A/cm2]');
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