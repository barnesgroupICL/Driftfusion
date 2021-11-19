function TPV_struct = TPV_script(structs, pulse_duration, pulse_int_fraction, save_solutions)
%TPVCONST_FULL_EXEC - Simulate Transient PhotoVoltage in a range of background light intensities without varying the pulse intensity
% The pulse intensity is tuned to be a small perturbation compared to the 1
% sun background illumination, then the TPV is simulated for all the
% provided solutions without changing the chosen pulse intensity.
%
% Syntax:  TPV_struct = TPVconst_full_exec(SYMSTRUCTS, SAVE_SOLUTIONS, SAVE_RESULTS)
%
% Inputs:
%   SYMSTRUCTS - can be a cell structure containing structs at various background
%     light intensities. This can be generated using genIntStructs.
%     Otherwise it can be a single struct as created by PINDRIFT.
%   SAVE_SOLUTIONS - logical, whether to assign in volatile base
%     workspace the calculated solutions of single TPV decays
%
% Outputs:
%   TPV_STRUCT - a struct containing the most important results of the simulation
%
% Example:
%   TPVconst_full_exec(genIntStructs(ssol_i_eq, 100, 0.1, 4, true), true, true)
%     calculate also with dark background
%   TPVconst_full_exec(ssol_i_light, true, true)
%     calculate TPV just for one given struct
%
% Other m-files required: TPV_single_exec, TPV_single_analysis, TPV_full_analysis
% Subfunctions: none
% MAT-files required: none
%
% See also genIntStructs, pindrift, TPV_single_exec, TPV_single_analysis.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

% in case a single struct is given in input, convert it to a cell structure
% with just one cell
if length(structs(:, 1)) == 1 % if the input is a single structure instead of a cell with structures
    symstructs_temp = cell(2, 1);
    symstructs_temp{1, 1} = structs;
    symstructs_temp{2, 1} = inputname(1);
    structs = symstructs_temp;
end

% don't display figures until the end of the script, as they steal the focus
% taken from https://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
set(0, 'DefaultFigureVisible', 'off');
min_tmax_to_pulse_duration_ratio = 1000;
tmax = max(1, pulse_duration*min_tmax_to_pulse_duration_ratio);
tpoints = 200;


%% pre allocate arrays filling them with zeros
Voc_array = zeros(length(structs(1, :)), 1);
delta_v_array = Voc_array;
monoexp_lin_array = Voc_array;
monoexp_tau_array = Voc_array;
biexp_lin_fast_array = Voc_array;
biexp_tau_fast_array = Voc_array;
biexp_lin_slow_array = Voc_array;
biexp_tau_slow_array = Voc_array;

%% do a serie of TPV measurements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This could be parallelized using parfor from Parallel Computing Toolbox %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(structs(1, :))
    bias_int = structs{1, i}.par.int1;
    % a pulse with 1 sun power would have pulse_int = 2.5, but in this case also pulse length has to be considered
    pulse_int = max(bias_int,1e-4) * pulse_int_fraction;
    disp([mfilename ' - Doing the TPV at ' num2str(bias_int) ' light intensity']);

    % sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
    sol_ill = lightonRs(structs{1, i}, bias_int, -10, true, 1e6, 10);

    warning('off', 'Driftfusion:verifyStabilization');
    reachedStability = false;
    while ~reachedStability
        sol_TPV = doTPVonly(sol_ill, pulse_int, tmax, tpoints, 100*pulse_duration/tmax, 2); % do single TPV
        % the time_fraction should be larger than the
        % 1/min_tmax_to_pulse_duration_ratio
        reachedStability = verifyStabilization(sol_TPV.u, sol_TPV.t, 0.002);
        tmax = 10 * tmax;
    end
    warning('on', 'Driftfusion:verifyStabilization');
    
    if save_solutions
        sol_name = matlab.lang.makeValidName([structs{2, i} '_TPV_' num2str(pulse_int_fraction) '_' num2str(pulse_duration)]);
        assignin('base', sol_name, sol_TPV);
    end
    tmax = max(sol_TPV.par.tmax/10, pulse_duration*10); % save the last good tmax for next cycle
    % extract parameters and do plot
    [Voc_array(i), delta_v_array(i), monoexp_lin_array(i),...
        monoexp_tau_array(i), biexp_lin_fast_array(i),...
        biexp_tau_fast_array(i), biexp_lin_slow_array(i),...
        biexp_tau_slow_array(i)] = TPV_ana(sol_TPV, false);
end

% the pendent of the high illumination region could be useful to calculate
% weight = (Voc_array-min(Voc_array))/(max(Voc_array)-min(Voc_array));
% F = @(coef,xdata)coef(1)*exp(xdata*coef(2));
% monoexp_pendent_fit = fitnlm(Voc_array, monoexp_tau_array, 'Weigth', weigth);
% pendent = monoexp_pendent_fit.Coefficients.Estimate(2);

%% save results
TPV_struct.sol_name = structs{2, 1};
TPV_struct.Voc = Voc_array;
TPV_struct.monoexp_tau = monoexp_tau_array;
TPV_struct.monoexp_lin = monoexp_lin_array;
TPV_struct.biexp_tau_fast = biexp_tau_fast_array;
TPV_struct.biexp_tau_slow = biexp_tau_slow_array;
TPV_struct.biexp_lin_fast = biexp_lin_fast_array;
TPV_struct.biexp_lin_slow = biexp_lin_slow_array;

%% plot results
TPV_script_ana(TPV_struct);
    
% make the figures appear, all at the end of the script
set(0, 'DefaultFigureVisible', 'on');
figHandles = findall(groot, 'Type', 'figure');
set(figHandles(:), 'visible', 'on')

%------------- END OF CODE --------------
