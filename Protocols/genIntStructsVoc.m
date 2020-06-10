function [goodVocAsymStructCell, VOCs] = genIntStructsVoc(struct_eq, startInt, endInt, points, include_dark)
% GENINTSTRUCTSREALVOC - Generates a cell containing structures of solutions at various light intensities at an accurate VOC
% This script just uses other three scripts: genIntStructs, findOptimVoc
% and asymmetricize. Both symmetric and asymmetric solutions are supported
% in input, but the usage of symmetric solutions is strongly encouraged as
% the findOptimVoc will start from a condition closer to the Voc.
%
% Syntax:  [structCell, VOCs] = genIntStructs(struct_eq, startInt, endInt, points, include_dark)
%
% Inputs:
%   STRUCT_EQ - a solution struct as created by DRIFTFUSION in dark
%     conditions, preferably a symmetric (open circuit) solution
%   STARTINT - higher requested illumination.
%   ENDINT - lower requested illumination.
%   POINTS - number of illumination requested between STARTINT and ENDINT, including extrema, except dark.
%   INCLUDE_DARK - logical, if to include the dark solution in the output
%     structure
%
% Outputs:
%   GOODVOCASYMSTRUCTCELL - a cell containing structs of asymmetric solutions at various light
%     intensities with an applied voltage equal to the VOC
%   VOCS - an array with the VOC, getting populated just if the input structures were at open circuit
%
% Example:
%   [structs_voc, VOCs] = genIntStructsRealVoc(ssol_i_eq_SR, 1, 1e-3, 7, true)
%     prepare solutions at 100, 10, 1, 0.1 and 0 illumination intensities
%
% Other m-files required: changeLight, pindrift, genIntStructs,
%   asymmetricize, findOptimVoc
% Subfunctions: none
% MAT-files required: none
%
% See also genVappStructs, changeLight, pindrift, genIntStructs, findOptimVoc.
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% May 2018; Last revision: May 2018

%------------- BEGIN CODE --------------

if ~struct_eq.par.OC
    disp([mfilename ' - the input solution is in short circuit! ' mfilename ' runs much faster when starting from OC conditions!']);
end

% use the normal genIntStructs for obtaining a cell with structs at various
% light intensities
[badVocStructCell, ~, ~] = genIntStructs(struct_eq, startInt, endInt, points, include_dark);

% how many solutions have been obtained
nsolutions = length(badVocStructCell(1,:));

% preallocate
VOCs = NaN(1, nsolutions);
% the second row of this cell, containing the solutions names, will not be
% changed
goodVocAsymStructCell = badVocStructCell;

%% generate solutions
for i = 1:nsolutions
    
    % decrease annoiance by figures popping up
    badVocStructCell{1, i}.par.figson = 0;
    
    asymstruct_Int = badVocStructCell{1, i};
    % the asymmetricized solution could require some stabilization after
    % breaking
    asymstruct_Int = stabilize(asymstruct_Int);
    
    % use findOptimVoc for finding the applied voltage that minimizes the
    % residual current
    disp([mfilename ' - finding real Voc for illumination intensity ' num2str(badVocStructCell{1, i}.par.int1)])
    light_intensity = asymstruct_Int.par.int1;
    [~, VOC_Rs] = findVocDirect(asymstruct_Int, light_intensity, 1);
    [asymstruct_Int_Voc, VOC] = findVoc(asymstruct_Int, light_intensity, asymstruct_Int.par.mobseti, VOC_Rs*0.95, VOC_Rs*1.05, 1e-5);

    % restore figson before saving
    asymstruct_Int_Voc.par.figson = 1;
    % replace the solution at the bad VOC with the new one
    goodVocAsymStructCell{1, i} = asymstruct_Int_Voc;
    % populate the array containing VOCs
    VOCs(i) = VOC;
end

%------------- END OF CODE --------------
