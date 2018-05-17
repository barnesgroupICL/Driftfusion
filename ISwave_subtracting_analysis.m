function [subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t] = ISwave_subtracting_analysis(asymstruct_ISwave)
%ISWAVE_SUBTRACTING_ANALYSIS - Calculate the charge excess in a device
% under illumination compared with the previous time point and express this
% as a current, for Impedance Spectroscopy with oscillating voltage (ISwave)
% This is used as a reference for IS and for separating the ionic contribution
%
% Syntax:  [subtracting_n_t, subtracting_n_intr_t, subtracting_n_contacts_t, subtracting_i_abs_t, subtracting_i_t] = ISwave_subtracting_analysis(asymstruct_ISwave)
%
% Inputs:
%   ASYMSTRUCT_ISWAVE - a struct with a solution being perturbated by an
%     oscillating voltage, as generated from ISwave_EA_single_exec
%
% Outputs:
%
%
% Example:
%   ISwave_subtracting_analysis(ISwave_EA_single_exec(asymmetricize(ssol_i_light), 2e-3, 1e6, 20, 40, true, true, 1e-4))
%     extract reference values
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also ISwave_full_exec, ISwave_EA_single_exec, ISwave_single_analysis.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

% evil shortcut
s = asymstruct_ISwave;

% Intrinsic points logical array
itype_points= (s.x >= s.p.tp & s.x <= s.p.tp + s.p.ti);

% make a new matrix subtracting from each profile at a time point the
% profile at the previous time point, the resulting rows will be one row
% less than the original solution matrix

% now I'm doing the gradient and then integrating it in the other
% dimension, it would be more efficient to do the integration before
% differentiating
[~, n_delta_matrix] = gradient(s.sol(:, :, 1), s.x, s.t);%s.sol(2:end, :, 1) - s.sol(1:end-1, :, 1);
[~, i_delta_matrix] = gradient(s.sol(:, :, 3), s.x, s.t);%s.sol(2:end, :, 3) - s.sol(1:end-1, :, 3);
i_delta_abs_matrix = abs(i_delta_matrix); % out of the intrinsic the difference is constantly zero

% sum of the electrons over thickness, obtain a point for each time
n_delta_t = transpose(trapz(s.x, n_delta_matrix, 2));

% get the middle of the intrinsic material, for this calculation maybe the
% middle could be not the correct point to use, but good enough
% get the index of the point closest to thickness p-type + half of
% thickness intrinsic
[~, x_center_index] = min(abs(s.x - (s.p.tp + s.p.ti / 2)));
% invert in sign the second spatial half of the ion profile, so that the
% integral over the thickness doesn't result always zero
i_delta_reflex_matrix = i_delta_matrix;
i_delta_reflex_matrix(:, x_center_index:end) = -i_delta_reflex_matrix(:, x_center_index:end);
% dividing by two for not counting twice each ion (once from where it left and once for where it arrived)
i_delta_t = transpose(trapz(s.x, i_delta_reflex_matrix, 2) / 2);

% sum of the absolute delta ions profile over thickness, obtain a point
% for each time, this is equivalent of a total current due to ions
% displacing
% dividing by two for not counting twice each ion (once from where it left and once for where it arrived)
i_delta_abs_t = transpose(trapz(s.x, i_delta_abs_matrix, 2) / 2);

% sum of the electrons just in the intrinsic layer, obtain a point for each time
n_delta_intr_t = transpose(trapz(s.x(itype_points), n_delta_matrix(:, itype_points), 2));

%% dividing by time delta converts the output in a current
% convert from number of electrons to charge per time point
subtracting_n_t = -s.p.e * n_delta_t;
% convert from number of ions to charge per time point
subtracting_i_abs_t = s.p.e * i_delta_abs_t;
% convert from number of ions to charge per time point
subtracting_i_t = s.p.e * i_delta_t;
% convert from number of electrons in the intrinsic to charge per time point
subtracting_n_intr_t = -s.p.e * n_delta_intr_t;
% this is the same as actually doing the integral in the contacts
subtracting_n_contacts_t = subtracting_n_t - subtracting_n_intr_t;

%------------- END OF CODE --------------
