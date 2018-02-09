function structCell = genVappStructs(asymstruct, startVapp, endVapp, points)
%GENVAPPSTRUCTS - Generates a cell containing asymmetric structures of solutions at various applied voltages
%
% Syntax:  structCell = genVappStructs(asymstruct, startVapp, endVapp, points)
%
% Inputs:
%   ASYMSTRUCT - a solution asymmetric struct as created by PINDRIFT.
%   STARTVAPP - higher requested Vapp.
%   ENDVAPP - lower requested Vapp.
%   POINTS - number of applied voltages requested between STARTVAPP and ENDVAPP, including extrema.
%
% Outputs:
%   STRUCTCELL - a cell containing structs of solutions at various applied
%     voltages
%
% Example:
%   structs_Vapp_dark = genVappStructs(sol_i_eq, 1, 0, 6);
%     generates dark solutions at 6 different applied voltages
%
% Other m-files required: pindrift
% Subfunctions: none
% MAT-files required: none
%
% See also genIntStructs, pindrift.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

asymstruct_Vapp = asymstruct;

% estimate a good tmax
if asymstruct_Vapp.params.mui
    tmax_temp = 2^(-log10(asymstruct_Vapp.params.mui)) / 10 + 2^(-log10(asymstruct_Vapp.params.mue_i));
else
    tmax_temp = 2^(-log10(asymstruct_Vapp.params.mue_i));
end

% define linearly spaced applied voltage values
Vapp_array = linspace(startVapp, endVapp, points);
% usually the solution in short circuit is given in input, so start from
% zero Vapp
Vapp_array = sort(Vapp_array, 'ascend');

% pre allocate
structCell = cell(2, length(Vapp_array));

existingVars = evalin('base', 'who');

%% generate solutions
for i = 1:length(Vapp_array)
    disp([mfilename ' - applied voltage ' num2str(Vapp_array(i))])
    name = matlab.lang.makeValidName([inputname(1) '_Vapp_' num2str(Vapp_array(i))]);
    if any(strcmp(existingVars, name)) % check if a structure with the same name already exists
        asymstruct_Vapp = evalin('base', name);
        % if none of the previous if conditions is true, just use the
        % previously set struct_Vapp
    end
    % decrease annoiance by figures popping up
    asymstruct_Vapp.params.figson = 0;
    
    p = asymstruct_Vapp.params;
    % prepare parameters for the voltage change
    p.tmesh_type = 2;
    p.t0 = 1e-10;
    p.tpoints = 50; % wide, does not matter much but have to be bigger than t_npoints_V_step
    p.calcJ = 0; % we do not need to calculate current here
    p.tmax = tmax_temp;

    % ideal sudden change is voltage maybe cannot be solved
    % boundary conditions have to be changed more slowly 
    % the solver needs at least one point at the starting voltage, before switching to short circuit
    p.JV = 2; % new mode for fast Voc changes
    p.t_npoints_V_step = 10; % new variable for Voc sudden changes, if equals tpoints is like a JV scan, minimum 2 points which is a direct step
    p.Vstart = asymstruct_Vapp.Efn(end) - asymstruct_Vapp.Efp(1);
    p.Vend = Vapp_array(i); % final Voc

    % go to the new Vapp
    asymstruct_Vapp = pindrift(asymstruct_Vapp, p);

    % restore figson before saving
    asymstruct_Vapp.params.figson = 1;
    structCell{1, i} = asymstruct_Vapp;
    structCell{2, i} = name;
    assignin('base', name, asymstruct_Vapp);
end

%------------- END OF CODE --------------
