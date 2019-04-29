function [structCell, V_array, J_array] = genIntStructs(struct_eq, startInt, endInt, points, include_dark)
%GENINTSTRUCTS - Generates a cell containing structures of solutions at various light intensities, starting from the one in dark
%
% Syntax:  [structCell, V_array, J_array] = genIntStructs(struct_eq, startInt, endInt, points, include_dark)
%
% Inputs:
%   STRUCT_EQ - a solution struct as created by df in dark conditions, or a logic false for not having the dark.
%   STARTINT - higher requested illumination.
%   ENDINT - lower requested illumination.
%   POINTS - number of illumination requested between STARTINT and ENDINT, including extrema, except dark.
%   INCLUDE_DARK - logical, if to include the dark solution in the output
%     structure
%
% Outputs:
%   STRUCTCELL - a cell containing structs of solutions at various light
%     intensities
%   V_ARRAY - an array with the voltages present in the solution, which is
%     an aproximation of the VOCs, getting populated just if the input
%     structures were at open circuit
%   J_ARRAY - an array with the currents present in the solutionJs
%
% Example:
%   [structs_oc, VOCs, ~] = genIntStructs(ssol_i_eq_SR, 10, 0.01, 4, true)
%     prepare open circuit solutions at 10, 1, 0.1, 0.01 and dark illumination intensities
%   [structs_sc, ~, JSCs] = genIntStructs(sol_i_eq_SR, 100, 0.1, 4, false)
%     prepare short circuit solutions at 100, 10, 1 and 0.1 illumination intensities
%
% Other m-files required: changeLight, df, dfana
% Subfunctions: none
% MAT-files required: none
%
% See also genVappStructs, changeLight, df.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

struct_Int = struct_eq;

% define light intensity values, always decreasing
Int_array = logspace(log10(startInt), log10(endInt), points);
% if include_dark, the input dark solution is included as-it-comes in the
% output structure
if include_dark
    Int_array = [Int_array, 0];
end
% decrease annoiance by figures popping up
struct_eq.par.figson = 0;

Int_array = sort(Int_array, 'ascend'); % from high intensity to low, beware to not use changeLight for reaching dark as logspace to zero doesn't work

% pre allocate
structCell = cell(2, length(Int_array));
V_array = NaN(1, length(Int_array));
J_array = NaN(1, length(Int_array));

changeLight_tmax = false; % automatic tmax for the first run

%% generate solutions
for i = 1:length(Int_array)
    disp([mfilename ' - illumination intensity ' num2str(Int_array(i))])
    name = matlab.lang.makeValidName([inputname(1) '_Int_' num2str(Int_array(i))]);
    % decrease annoiance by figures popping up
    struct_Int.par.figson = 0;
    % in any case except dark (it would work, should not be needed), stabilize at the new intensity
    if Int_array(i)
        struct_Int = changeLight(struct_Int, Int_array(i), changeLight_tmax); % change light intensity
        changeLight_tmax = max(struct_Int.par.tmax / 5, 1e-8); % time to use for next iteration
    end

    [V_temp, J_temp] = dfana(struct_Int);
    V_array(i) = V_temp(end);
    J_array(i) = J_temp(end);

    % restore figson before saving
    struct_Int.par.figson = 1;
    structCell{1, i} = struct_Int;
    structCell{2, i} = name;
end

%------------- END OF CODE --------------
