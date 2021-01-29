function [structs_oc, VOCs] = genIntStructsVoc(struct_eq, startInt, endInt, points, include_dark)
% GENINTSTRUCTSVOC - Generates a cell containing structures of solutions at various light intensities at an accurate VOC
% This script just uses other two scripts: genIntStructs, findVocOptim.
%
% Syntax:  [structs_oc, VOCs] = genIntStructsVoc(struct_eq, startInt, endInt, points, include_dark)
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
%   STRUCTS_OC - a cell containing structs of asymmetric solutions at various light
%     intensities with an applied voltage equal to the VOC
%   VOCS - an array with the calculated VOC
%
% Example:
%   [structs_oc, VOCs] = genIntStructsVoc(soleq.ion, 10, 1e-3, 5, true)
%     prepare solutions at 10, 1, 0.1, 0.01, 0.001 and 0 illumination intensities
%
% Other m-files required: df, genIntStructs, findVocOptim
% Subfunctions: none
% MAT-files required: none
%
% See also genVappStructs, changeLight, df, genIntStructs, findVocOptim.
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

% use the normal genIntStructs for obtaining a cell with structs at various
% light intensities
[structs_sc, ~, ~] = genIntStructs(struct_eq, startInt, endInt, points, include_dark);

% how many solutions have been obtained
nsolutions = length(structs_sc(1,:));

% preallocate
VOCs = NaN(1, nsolutions);
% the second row of this cell, containing the solutions names, will not be
% changed
structs_oc = structs_sc;

%% generate solutions
for i = 1:nsolutions
    
    % decrease annoiance by figures popping up
    %SCStructCell{1, i}.par.figson = 0;
    
    struct_Int = structs_sc{1, i};
    % the asymmetricized solution could require some stabilization after
    % breaking
    %struct_Int = stabilize(struct_Int);
    
    % use findVocOptim for finding the applied voltage that minimizes the
    % residual current
    disp([mfilename ' - finding real Voc for illumination intensity ' num2str(structs_sc{1, i}.par.int1)])

    % if solution seems at short circuit, estimate a voltage based on int1
    % illumination
    if abs(struct_Int.par.Vapp) < 0.1
        % findVocDirect could be used but it's slow
        %[~, guessVoc] = findVocDirect(struct_Int, struct_Int.par.int1, 1);

        % very rough estimation of guessVoc
        guessVoc = 1 + 0.035 * log(struct_Int.par.int1 - 0.0002);
        [struct_Int_Voc, VOC] = findVocOptim(struct_Int, guessVoc);
    else
        [struct_Int_Voc, VOC] = findVocOptim(struct_Int);
    end

    % replace the solution at the bad VOC with the new one
    structs_oc{1, i} = struct_Int_Voc;
    % populate the array containing VOCs
    VOCs(i) = VOC;
end

%------------- END OF CODE --------------
