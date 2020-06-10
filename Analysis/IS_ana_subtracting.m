function [subtracting_q_t, subtracting_q_intr_t] = IS_ana_subtracting(struct_IS)
%IS_ANA_SUBTRACTING - Calculates the time derivative of the total charge in the device
% Calculate the charge excess in a device compared with the previous time point and express this as a current
% for Impedance Spectroscopy with oscillating voltage (IS)
% This is used as a reference for IS and for separating the ionic contribution
%
% Syntax:  [subtracting_n_t, subtracting_n_intr_t] = IS_ana_subtracting(struct_IS)
%
% Inputs:
%   STRUCT_IS - a struct with a solution being perturbed by an
%     oscillating voltage, as generated from doIS_EA
%
% Outputs:
%   subtracting_n_t - charge variation over time in the whole device
%   subtracting_n_intr_t - charge variation over time in the perovskite
%     layer
%
% Example:
%   [subtracting_n_t, subtracting_n_intr_t] = IS_ana_subtracting(ssol_i_1S_SR_is_100mHz_2mV)
%     extract reference values
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also IS_script, doIS_EA, IS_ana.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: May 2020

%------------- BEGIN CODE --------------

% shortcut
s = struct_IS;

% Intrinsic points logical array
itype_points = logical(s.par.dev.g0);

% make a new matrix subtracting from each profile at a time point the
% profile at the previous time point, the resulting rows will be one row
% less than the original solution matrix

% now I'm doing the gradient and then integrating it in the other
% dimension, it would be more efficient to do the integration before
% differentiating
[~, q_delta_matrix] = gradient((s.u(:, :, 2)+s.u(:, :, 3))/2, s.x, s.t);

% sum of the electrons over thickness, obtain a point for each time
q_delta_t = transpose(trapz(s.x, q_delta_matrix, 2));

% sum of the electrons just in the intrinsic layer, obtain a point for each time
n_delta_intr_t = transpose(trapz(s.x(itype_points), q_delta_matrix(:, itype_points), 2));

%% dividing by time delta converts the output in a current
% convert from number of electrons to charge per time point
subtracting_q_t = s.par.e * q_delta_t;
% convert from number of electrons in the intrinsic to charge per time point
subtracting_q_intr_t = s.par.e * n_delta_intr_t;

%------------- END OF CODE --------------
