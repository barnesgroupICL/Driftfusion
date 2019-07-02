function CE_results = CE_full_exec(symstructs, do_graphics, save_solutions)
%CE_FULL_EXEC - Do Charge Extraction (CE) in a range of background light intensities
%
% Syntax:  CE_results = CE_full_exec(symstructs, save_solutions, save_results)
%
% Inputs:
%   SYMSTRUCTS - can be a cell structure containing structs at various background
%     light intensities. This can be generated using genIntStructs.
%     Otherwise it can be a single struct as created by PINDRIFT.
%   DO_GRAPHICS - logical, whether to graph the individual solutions and
%     the overall graphics
%   SAVE_SOLUTIONS - logical, whether to assign in volatile base
%     workspace the calculated solutions of single CE decays
%
% Outputs:
%   CE_results - a struct containing the most important results of the simulation
%
% Example:
%   CE_full_exec(genIntStructs(ssol_i_eq_SR, 100, 1e-7, 4, false), true, true)
%     calculate also with dark background
%   CE_full_exec(ssol_i_light, true, true)
%     calculate CE just for one given struct
%
% Other m-files required: CE_single_exec, CE_ISstep_single_analysis,
%   CE_ISstep_subtracting_analysis, CE_full_fit, CE_full_analysis,
%   asymmetricize, pinana
% Subfunctions: none
% MAT-files required: none
%
% See also genIntStructs, pindrift, CE_single_exec, CE_ISstep_single_analysis, CE_ISstep_subtracting_analysis, CE_full_fit, CE_full_analysis.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

% increase graphics font size
set(0, 'defaultAxesFontSize', 20)

% in case a single struct is given in input, convert it to a cell structure
% with just one cell
if length(symstructs(:, 1)) == 1 % if the input is a single structure instead of a cell with structures
    symstructs_temp = cell(2, 1);
    symstructs_temp{1, 1} = symstructs;
    symstructs_temp{2, 1} = inputname(1);
    symstructs = symstructs_temp;
end

% don't display figures until the end of the script, as they steal the focus
% taken from https://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
set(0, 'DefaultFigureVisible', 'off');

tmax = 1e-3; % small tmax for starting

%% pre allocate arrays filling them with zeros
Voc_array = zeros(length(symstructs(1, :)), 1);
Int_array = Voc_array;
extracting_charges_array = Voc_array;
early_extracting_charges_array = Voc_array;
subtracting_charges_array = Voc_array;
subtracting_charges_contacts_array = Voc_array;
subtracting_charges_intrinsic_array = Voc_array;
ions_displaced_array = Voc_array;

%% do a serie of CE measurements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This could be parallelized using parfor from Parallel Computing Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([mfilename ' - Doing the CE at various light intensities']);
for i = 1:length(symstructs(1, :))
    Int_array(i) = symstructs{1, i}.p.Int;
    % decrease annoiance by figures popping up
    symstructs{1, i}.p.figson = 0;
    asymstruct_Int = asymmetricize(symstructs{1, i});
    % the asymmetricized solution often (mostly when both high illumination
    % and mobile ions are present) needs stabilization
    asymstruct_Int = stabilize(asymstruct_Int);
    [Vapp_arr, ~, ~] = pinana(asymstruct_Int);
    Voc_array(i) = Vapp_arr(end);
    asymstruct_Int_CE = CE_single_exec(asymstruct_Int, tmax); % do CE
    if save_solutions
        sol_name = matlab.lang.makeValidName([symstructs{2, i} '_CE']);
        asymstruct_Int_CE.p.figson = 1; % re-enable figures
        assignin('base', sol_name, asymstruct_Int_CE);
    end
    % tmax can be a bit shorter with lower light intensities, doesn't
    % really change much, anyway check if can be halved
    warning('off', 'pindrift:verifyStabilization'); % here the stabilization is not expected nor needed, disable warning
    if verifyStabilization(asymstruct_Int_CE.sol, asymstruct_Int_CE.t, 0.1)
        tmax = asymstruct_Int_CE.p.tmax / 2;
    else
        tmax = asymstruct_Int_CE.p.tmax;
    end
    warning('on', 'pindrift:verifyStabilization');
    
    % extract parameters and do plot
    [temp_charges_array, index_early] = CE_ISstep_single_analysis(asymstruct_Int_CE, Int_array(i), false, do_graphics);
    extracting_charges_array(i) = temp_charges_array(end); % last value of charge integration
    early_extracting_charges_array(i) = temp_charges_array(index_early); % early value of charge integration
    [subtracting_charges_array(i), subtracting_charges_intrinsic_array(i)] =...
        CE_ISstep_subtracting_analysis(asymstruct_Int, asymstruct_Int_CE); % extract parameters via solution comparison and do plot
end

%% fit the charge vs voltage points with a linear + exponential curve
extracting_cap_array = CE_full_fit(Voc_array, extracting_charges_array);
early_cap_array = CE_full_fit(Voc_array, early_extracting_charges_array);
subtracting_charges_contacts_array = subtracting_charges_array - subtracting_charges_intrinsic_array;
subtracting_charges_intrinsic_cap_array = CE_full_fit(Voc_array, subtracting_charges_intrinsic_array);

%% save results
sun_index = find(Int_array == 1); % used for plotting vertical lines at 1 sun

CE_results.sol_name = symstructs{2, 1};
CE_results.Voc = Voc_array;
CE_results.sun_index = sun_index;
CE_results.extracting_charges = extracting_charges_array;
CE_results.early_extracting_charges = early_extracting_charges_array;
CE_results.subtracting_charges = subtracting_charges_array;
CE_results.subtracting_charges_contacts = subtracting_charges_contacts_array;
CE_results.subtracting_charges_intrinsic = subtracting_charges_intrinsic_array;
%CE_results.ions_displaced = ions_displaced_array;
CE_results.extracting_geom_cap = extracting_cap_array(1);
CE_results.extracting_chem_cap_amp = extracting_cap_array(2);
CE_results.extracting_chem_cap_gamma = extracting_cap_array(3);
CE_results.early_geom_cap = early_cap_array(1);
CE_results.early_chem_cap_amp = early_cap_array(2);
CE_results.early_chem_cap_gamma = early_cap_array(3);
CE_results.subtracting_intrinsic_geom_cap = subtracting_charges_intrinsic_cap_array(1);
CE_results.subtracting_intrinsic_chem_cap_amp = subtracting_charges_intrinsic_cap_array(2);
CE_results.subtracting_intrinsic_chem_cap_gamma = subtracting_charges_intrinsic_cap_array(3);
CE_results.BC = symstructs{1, 1}.p.BC;
CE_results.calcJ = symstructs{1, 1}.p.calcJ;


%% plot results
if do_graphics
    CE_full_analysis(CE_results);
end

% make the figures appear, all at the end of the script
set(0, 'DefaultFigureVisible', 'on');
figHandles = findall(groot, 'Type', 'figure');
set(figHandles(:), 'visible', 'on')

%------------- END OF CODE --------------
