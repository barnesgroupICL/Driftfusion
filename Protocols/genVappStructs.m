function structCell = genVappStructs(asymstruct, Vapp_array)
%GENVAPPSTRUCTS - Generates a cell containing asymmetric structures of solutions at various applied voltages
%
% Syntax:  structCell = genVappStructs(asymstruct, Vapp_array)
%
% Inputs:
%   ASYMSTRUCT - a solution asymmetric struct as created by PINDRIFT.
%   VAPP_ARRAY - an array containing the requested Vapp list.
%
% Outputs:
%   STRUCTCELL - a cell containing structs of solutions at various applied
%     voltages, ordered with ascending voltages
%
% Example:
%   structs_Vapp_dark = genVappStructs(sol_i_eq_SR, 0:0.2:1);
%     generates dark solutions at 0, 0.2, 0.4, 0.6, 0.8 and 1 V applied voltages
%   structs_Vapp_dark = genVappStructs(sol_i_eq_SR, linspace(0, 1, 15));
%     generates dark solutions at 15 different voltages from 0 to 1 V
%   structs_Vapp_dark = genVappStructs(sol_i_eq_SR, [0, 0.8, 0.9]);
%     generates dark solutions at a list of applied voltages, remember that
%     the list will have an ascending ordering in the output structure
%
% Other m-files required: df, dfana
% Subfunctions: none
% MAT-files required: none
%
% See also genIntStructs, df.

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
if asymstruct_Vapp.par.mucat
    tmax_temp = min(1e1, 2^(-log10(asymstruct_Vapp.par.mucat)) / 10 + 2^(-log10(asymstruct_Vapp.par.mue(1))));
else
    tmax_temp = min(1e-2, 2^(-log10(asymstruct_Vapp.par.mue(1))));
end

% usually the solution in short circuit is given in input, so start from
% zero Vapp
Vapp_array = sort(Vapp_array, 'ascend');

% pre allocate
structCell = cell(2, length(Vapp_array));

%% generate solutions
for i = 1:length(Vapp_array)
    disp([mfilename ' - applied voltage ' num2str(Vapp_array(i))])
    name = matlab.lang.makeValidName([inputname(1) '_Vapp_' num2str(Vapp_array(i))]);
    % decrease annoiance by figures popping up
    asymstruct_Vapp.par.figson = 0;
    
    par = asymstruct_Vapp.par;
    % prepare parameters for the voltage change
    par.tmesh_type = 2;
    par.t0 = 1e-10;
    par.tpoints = 50; % wide, does not matter much but have to be bigger than t_npoints_V_step
    par.calcJ = 0; % we do not need to calculate current here
    par.tmax = tmax_temp;

    % ideal sudden change is voltage maybe cannot be solved
    % boundary conditions have to be changed more slowly 
    % the solver needs at least one point at the starting voltage, before switching to short circuit
    par.JV = 2; % new mode for arbitrary Vapp profiles
    [Vapp_arr, ~] = dfana(asymstruct_Vapp); % quasi-fermi levels are needed for knowing the current Voc
    Vstart = Vapp_arr(end); % current Voc
    Vend = Vapp_array(i); % final Voc
    
    % the third parameter is the time point when the change from the old
    % voltage to the new one happens
    par.Vapp_params = [Vstart, Vend, 10 * par.t0];
    % starts at Vstart, linear Vapp decrease for the first points, then constant to Vend
    par.Vapp_func = @(coeff, t) (coeff(1) + (coeff(2)-coeff(1)).*t/coeff(3)) .* (t < coeff(3)) + coeff(2) .* (t >= coeff(3));

    % go to the new Vapp
    asymstruct_Vapp = df(asymstruct_Vapp, par);
    
    %% stabilize
    
    % should be set anyway, but for avoiding risks, set Vapp
    asymstruct_Vapp.par.Vapp = Vend;
    % eliminate JV configuration before saving (useless as stabilize would eliminate it anyway...)
    asymstruct_Vapp.par.JV = 0;
    
    asymstruct_Vapp = stabilize(asymstruct_Vapp); % go to steady state
    
    % restore figson before saving
    asymstruct_Vapp.par.figson = 1;

    % re-establish the original calcJ
    asymstruct_Vapp.par.calcJ = asymstruct.par.calcJ;

    structCell{1, i} = asymstruct_Vapp;
    structCell{2, i} = name;
end

%------------- END OF CODE --------------
