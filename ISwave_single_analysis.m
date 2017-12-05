function [fit_coeff, fit_idrift_coeff, subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t, Ji_disp] = ISwave_single_analysis(sol_i_Int_ISwave)
% calculate impedance (reactance and resistance) and phase by impedance
% spectroscopy with oscillating voltage

% evil shortcut
s = sol_i_Int_ISwave;

[subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t, Ji_disp] = ISwave_subtracting_analysis(sol_i_Int_ISwave);

Vapp = s.params.Vapp_func(s.params.Vapp_params, s.t);

% round should not be needed here
% s.params.Vapp_params(4) should be pulsatance
periods = round(s.params.tmax * s.params.Vapp_params(4) / (2 * pi));
% this works if, as now is, the frequency is constant internally each solution
% round should not be needed here
fit_t_index = round((s.params.tpoints - 1) * round(periods * 0.5) / periods);
fit_t = s.t(fit_t_index:end);
fit_J = s.Jtotr(fit_t_index:end) / 1000; % in Ampere
fit_Jidrift = s.Jidrift_points(fit_t_index:end, end) / 1000; % in Ampere
fit_coeff = ISwave_single_fit(fit_t, fit_J, s.params.J_func);
if isfield(s.params, 'Jpoints') % if the ionic drift has been calculated, fit it
    fit_idrift_coeff = ISwave_single_fit(fit_t, fit_Jidrift, s.params.J_func);
end

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
%     plot(s.t, s.Jpartr);
    plot(s.t, s.Jdispr / 1000); % Ampere
%     plot(s.t(2:end), subtracting_n_t);
    plot(s.t(2:end), subtracting_n_intr_t);
    plot(s.t(2:end), subtracting_n_contacts_t);
    plot(s.t(2:end), subtracting_i_abs_t);
    plot(s.t(2:end), subtracting_i_t);
    plot(s.t, Ji_disp);
    plot(fit_t, s.params.J_func(fit_coeff, fit_t), 'g')
    plot(fit_t, J_outphase, 'r')
    legend_array = ["Total J", "Displacement J", "Charges intrinsic", "Charges contacts", "Displaced ions abs", "Displaced ions", "Disp Curr Ions", "Fit", "Out of Phase Current"];
    if isfield(s.params, 'Jpoints')
%         plot(s.t, s.Jidiff_points(:,end), 'k');
        plot(s.t, s.Jidrift_points(:,end) / 1000, 'k.-'); % in Ampere
        plot(fit_t, s.params.J_func(fit_idrift_coeff, fit_t), 'g--');
        legend_array = [legend_array, "Ionic drift", "Ion drift fit"];
    end
    ylabel('Current [A/cm2]');
    hold off
    xlabel('Time [s]');
    legend(legend_array);
    hold off
