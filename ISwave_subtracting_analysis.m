function [subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_a_t] = ISwave_subtracting_analysis(sol_i_Int_ISwave)
% calculate the charge excess in device under illumination compared with
% the previous time point and express this as a current, for impedance
% spectroscopy with oscillating voltage, this is used as a reference for
% Impedance Spectroscopy and for separating the ionic contribution

% evil shortcut
s = sol_i_Int_ISwave;

% make a new matrix subtracting from each profile at a time point the
% profile at the previous time point, the resulting rows will be one row
% less than the original solution matrix

% for simulating a current, taking half of the abs is the right way
% I divide by two 'cause otherwise I would calculate once the charge getting to the new place and another time the charge missing in the old place
n_delta_matrix = abs(s.sol(2:end, :, 1) - s.sol(1:end-1, :, 1)) / 2;
a_delta_matrix = abs(s.sol(2:end, :, 3) - s.sol(1:end-1, :, 3)) / 2; % out of the intrinsic the difference is constantly zero

% t mesh should be equally spaced for IS wave experiment, anyway better to
% use actual time points just in case this will change in the future
t_array = s.t(2:end) - s.t(1:end-1);

% sum of the electrons over thickness, obtain a point for each time
n_delta_t = transpose(trapz(s.x, n_delta_matrix, 2));
a_delta_t = transpose(trapz(s.x, a_delta_matrix, 2));
% p-type points logical array
%ptype_points = (s.x <= s.params.tp);
% Intrinsic points logical array
itype_points= (s.x > s.params.tp & s.x < s.params.tp + s.params.ti);
% n-type points logical array
%ntype_points = (s.x >= s.params.tp + s.params.ti & s.x <= s.params.tp + s.params.ti + s.params.tn);
n_delta_intr_t = transpose(trapz(s.x(itype_points), n_delta_matrix(:, itype_points), 2));
assignin('base', 'sxi', s.x(itype_points))
assignin('base', 'n_delta_matrixi', n_delta_matrix(:, itype_points))
assignin('base', 'sx', s.x)
assignin('base', 'n_delta_matrix', n_delta_matrix)


% dividing by time delta converts the output in a current
subtracting_n_t = -s.params.e * n_delta_t ./ t_array; % convert from number of electrons to charge per time point
subtracting_a_t = s.params.e * a_delta_t ./ t_array; % convert from number of ions to charge per time point
subtracting_n_intr_t = -s.params.e * n_delta_intr_t ./ t_array; % convert from number of electrons in the intrinsic to charge per time point
subtracting_n_contacts_t = subtracting_n_t - subtracting_n_intr_t; % this is the same as actually doing the integral in the contacts