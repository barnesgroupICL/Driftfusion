function [subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t, Ji_disp] = ISwave_subtracting_analysis(sol_i_Int_ISwave)
% calculate the charge excess in device under illumination compared with
% the previous time point and express this as a current, for impedance
% spectroscopy with oscillating voltage, this is used as a reference for
% Impedance Spectroscopy and for separating the ionic contribution

% evil shortcut
s = sol_i_Int_ISwave;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select x points in each layer

% p-type points logical array
%ptype_points = (s.x <= s.params.tp);
% Intrinsic points logical array
itype_points= (s.x > s.params.tp & s.x < s.params.tp + s.params.ti);
% n-type points logical array
%ntype_points = (s.x >= s.params.tp + s.params.ti & s.x <= s.params.tp + s.params.ti + s.params.tn);

% make a new matrix subtracting from each profile at a time point the
% profile at the previous time point, the resulting rows will be one row
% less than the original solution matrix

% for simulating a current, taking half of the abs is the right way
% I divide by two 'cause otherwise I would calculate once the charge getting to the new place and another time the charge missing in the old place
n_delta_matrix = s.sol(2:end, :, 1) - s.sol(1:end-1, :, 1);
i_delta_matrix = s.sol(2:end, :, 3) - s.sol(1:end-1, :, 3);
i_delta_abs_matrix = abs(s.sol(2:end, :, 3) - s.sol(1:end-1, :, 3)); % out of the intrinsic the difference is constantly zero

% t mesh should be equally spaced for IS wave experiment, anyway better to
% use actual time points just in case this will change in the future
t_array = s.t(2:end) - s.t(1:end-1);

% sum of the electrons over thickness, obtain a point for each time
n_delta_t = transpose(trapz(s.x, n_delta_matrix, 2));

% get the middle of the intrinsic material, for this calculation maybe the
% middle could be not the correct point to use, but good enough
% get the index of the point closest to thickness p-type + half of
% thickness intrinsic
[~, x_center_index] = min(abs(s.x - (s.params.tp + s.params.ti / 2)));
% invert in sign the second spatial half of the ion profile, so that the
% integral over the thickness doesn't result always zero
i_delta_reflex_matrix = i_delta_matrix;
i_delta_reflex_matrix(:, x_center_index:end) = -i_delta_reflex_matrix(:, x_center_index:end);
% dividing by two for not counting twice each ion (once from where it left and once for where it arrived)
i_delta_t = transpose(trapz(s.x, i_delta_reflex_matrix, 2) / 2);

% subtract background ion concentration, for having less noise in trapz
i_matrix = s.sol(:, :, 3) - s.params.NI;
% calculate electric field due to ions in the middle of the material (no reason for taking the middle)
% Efield_i_middle = s.params.e * trapz(s.x(1:x_center_index), i_matrix(:, 1:x_center_index), 2) / s.params.eppi;
Efield_i = s.params.e * cumtrapz(s.x, i_matrix, 2) / s.params.eppi;
Efield_i_mean = mean(Efield_i(:, itype_points), 2);

Ji_disp = s.params.eppi * gradient(Efield_i_mean, s.t); % in Amperes

% sum of the absolute delta ions profile over thickness, obtain a point
% for each time, this is equivalent of a total current due to ions
% displacing
% dividing by two for not counting twice each ion (once from where it left and once for where it arrived)
i_delta_abs_t = transpose(trapz(s.x, i_delta_abs_matrix, 2) / 2);


% sum of the electrons just in the intrinsic layer, obtain a point for each time
n_delta_intr_t = transpose(trapz(s.x(itype_points), n_delta_matrix(:, itype_points), 2));

% dividing by time delta converts the output in a current
subtracting_n_t = -s.params.e * n_delta_t ./ t_array; % convert from number of electrons to charge per time point
subtracting_i_abs_t = s.params.e * i_delta_abs_t ./ t_array; % convert from number of ions to charge per time point
subtracting_i_t = s.params.e * i_delta_t ./ t_array; % convert from number of ions to charge per time point
subtracting_n_intr_t = -s.params.e * n_delta_intr_t ./ t_array; % convert from number of electrons in the intrinsic to charge per time point
subtracting_n_contacts_t = subtracting_n_t - subtracting_n_intr_t; % this is the same as actually doing the integral in the contacts