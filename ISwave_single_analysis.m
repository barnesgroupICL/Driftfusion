function [fit_coeff, subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_a_t] = ISwave_single_analysis(sol_i_Int_ISwave)
% calculate impedance (reactance and resistance) and phase by impedance
% spectroscopy with oscillating voltage

% evil shortcut
s = sol_i_Int_ISwave;

[subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_a_t] = ISwave_subtracting_analysis(sol_i_Int_ISwave);

Vapp = s.params.Vapp_func(s.params.Vapp_params, s.t);

% round should not be needed here
% s.params.Vapp_params(4) should be pulsatance
periods = round(s.params.tmax * s.params.Vapp_params(4) / (2 * pi));
% this works if, as now is, the frequency is constant internally each solution
% round should not be needed here
fit_t_index = round((s.params.tpoints - 1) * (periods - 1) / periods);
fit_t = s.t(fit_t_index:end);
fit_J = s.Jtotr(fit_t_index:end);
fit_coeff = ISwave_single_fit(fit_t, fit_J, s.params.Vapp_func);

% fourth value of Vapp_params is pulsatance
figure('Name', ['Single ISwave, Int ' num2str(s.params.Int) ' Freq ' num2str(s.params.Vapp_params(4) / (2 * pi))], 'NumberTitle', 'off');
    yyaxis right
    hold off
    plot(s.t, Vapp);
    ylabel('Applied voltage [V]');
    
    yyaxis left
    hold off
    plot(s.t, s.Jtotr);
    hold on
%     plot(s.t, s.Jpartr);
    plot(s.t, s.Jdispr);
%     plot(s.t(2:end), subtracting_n_t);
    plot(s.t(2:end), subtracting_n_intr_t * 1e3); % mA
    plot(s.t(2:end), subtracting_n_contacts_t * 1e3); % mA
    plot(s.t(2:end), subtracting_a_t * 1e3); % mA
    plot(fit_t, s.params.Vapp_func(fit_coeff, fit_t))
    legend_array = ["Total J", "Displacement J", "Charges intrinsic", "Charges contacts", "Displaced ions", "Fit"];
    if isfield(s, 'Jidiff_points')
%         plot(s.t, s.Jidiff_points(:,end), 'k');
        plot(s.t, s.Jidrift_points(:,end), 'k.-');
        legend_array = [legend_array, "Ionic drift"];
    end
    ylabel('Current [mA/cm2]');
    hold off
    xlabel('Time [s]');
    legend(legend_array);
    hold off