function [subtracting_charges, subtracting_charges_intrinsic] = CE_ISstep_subtracting_analysis(asymstruct_Int, asymstruct_CE_ISstep)
%CE_ISSTEP_SUBTRACTING_ANALYSIS - Calculate the charge excess in device under illumination
% The carriers density profile in the illuminated steady state solution is subtracted from the profiles
% at all time points of the solution evolving during the charge extraction experiment.
% These profiles are integrated over the device thickness, obtaining a time resolved excess charge amount.
% This is used as a reference for Charge Extraction and Impedance Spectroscopy with stepped voltage
%
% Syntax:  [subtracting_charges, subtracting_charges_contacts, subtracting_charges_intrinsic, ions_displaced] = CE_ISstep_subtracting_analysis(asymstruct_Int, asymstruct_CE_ISstep)
%
% Inputs:
%   ASYMSTRUCT_INT - struct of the stable solution before the CE or ISstep
%     experiment, to be used as a reference
%   ASYMSTRUCT_CE_ISSTEP - struct of the solution perturbed by CE or ISstep
%
% Outputs:
%   SUBTRACTING_CHARGES - array of excess charges in the whole device versus time
%   SUBTRACTING_CHARGES_INTRINSIC - array of excess charges in the intrinsic layer versus time
%
% Example:
%   CE_ISstep_subtracting_analysis(asymmetricize(ssol_i_light, 1), CE_single_exec(asymmetricize(ssol_i_light, 1), 1e-3, 1))
%     compare the solution before and after a CE experiment
%   CE_ISstep_subtracting_analysis(asymmetricize(ssol_i_light, 1), ISstep_single_exec(asymmetricize(ssol_i_light, 1), 1e-2, 1, 1e-3, 200))
%     compare the solution before and after an ISstep experiment
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also CE_single_exec, CE_full_exec, ISstep_single_exec, ISstep_full_exec.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

n_delta_profile = asymstruct_CE_ISstep.sol(end, :, 1) - asymstruct_Int.sol(end, :, 1);

% sum of the electrons as they could be extracted at the
% respective electrode completely
n_delta = trapz(asymstruct_Int.x, n_delta_profile);
subtracting_charges = -asymstruct_Int.p.e * n_delta; % convert from number of electrons to charge

% evil shortcut
s = asymstruct_Int;
% Intrinsic points logical array
itype_points= (s.x >= s.p.tp & s.x <= s.p.tp + s.p.ti);

n_delta_intrinsic = trapz(asymstruct_Int.x(itype_points), n_delta_profile(itype_points));
subtracting_charges_intrinsic = -asymstruct_Int.p.e * n_delta_intrinsic; % convert from number of electrons to charge

disp([mfilename ' - total electrons delta ' num2str(subtracting_charges)...
    ' C/cm2; intrinsic electrons delta ' num2str(subtracting_charges_intrinsic)...
    ' C/cm2.']);

%------------- END OF CODE --------------
