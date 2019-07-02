function TPV_struct = TPVconst_full_exec(symstructs, save_solutions)
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
if length(symstructs(:, 1)) == 1 % if the input is a single structure instead of a cell with structures
    symstructs_temp = cell(2, 1);
    symstructs_temp{1, 1} = symstructs;
    symstructs_temp{2, 1} = inputname(1);
    symstructs = symstructs_temp;
end

% don't display figures until the end of the script, as they steal the focus
% taken from https://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
set(0, 'DefaultFigureVisible', 'off');

% small tmax for getting a real delta_v, but enough for including also the
% peak, beware that t0 and pulselen are hardcoded in TPV-single-exec
tmax = 1e-5;
pulse_int = 100; % a pulse with 1 sun power would have pulse_int = 2.5, but in this case also pulse length has to be considered
delta_v_threshold = 3e-3; % for being sure of being in small perturbation regime

%% find laser pulse intensity for staying in the low perturbation regime
disp([mfilename ' - Reducing the pulse intensity for entering the small perturbation regime']);
warning('off','pindrift:verifyStabilization');
delta_v_ok = false;
% decrease annoiance by figures popping up
symstructs{1, 1}.p.figson = 0;

while ~delta_v_ok
    symstruct_pulse = TPV_single_exec(symstructs{1, 1}, 0, pulse_int);
    delta_v = max(symstruct_pulse.Voc); % the Voc baseline SHOULD BE set to zero
    if abs(min(symstruct_pulse.Voc)) > 0.5 * max(symstruct_pulse.Voc)
        warning('Big negative component in TPV, pulse intensity calibration will have wrong results');
    end
    disp([mfilename ' - pulse_int: ' num2str(pulse_int) '; delta_v: ' num2str(delta_v)]);
    if delta_v < delta_v_threshold
        delta_v_ok = true;
    else
        pulse_int = pulse_int * 0.98 * (delta_v_threshold / delta_v); % decrease pulse intensity until when is a small perturbation
    end
end
warning('on','pindrift:verifyStabilization');

%% pre allocate arrays filling them with zeros
Voc_array = zeros(length(symstructs(1, :)), 1);
Int_array = Voc_array;
delta_v_array = Voc_array;
monoexp_lin_array = Voc_array;
monoexp_tau_array = Voc_array;
biexp_lin_fast_array = Voc_array;
biexp_tau_fast_array = Voc_array;
biexp_lin_slow_array = Voc_array;
biexp_tau_slow_array = Voc_array;

%% do a serie of TPV measurements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This could be parallelized using parfor from Parallel Computing Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(symstructs(1, :))
    Int_array(i) = symstructs{1, i}.p.Int;
    disp([mfilename ' - Doing the TPV at ' num2str(Int_array(i)) ' light intensity']);
    % decrease annoiance by figures popping up
    symstructs{1, i}.p.figson = 0;

    symstruct_pulse = TPV_single_exec(symstructs{1, i}, tmax, pulse_int); % do single TPV
    if save_solutions
        sol_name = matlab.lang.makeValidName([symstructs{2, i} '_constpulse_' num2str(pulse_int)]);
        symstruct_pulse.p.figson = 1;
        assignin('base', sol_name, symstruct_pulse);
    end
    tmax = symstruct_pulse.p.tmax; % save the last good tmax for next cycle
    % extract parameters and do plot
    [Voc_array(i), delta_v_array(i), monoexp_lin_array(i),...
        monoexp_tau_array(i), biexp_lin_fast_array(i),...
        biexp_tau_fast_array(i), biexp_lin_slow_array(i),...
        biexp_tau_slow_array(i)] = TPV_single_analysis(symstruct_pulse);
end

% the pendent of the high illumination region could be useful to calculate
% weight = (Voc_array-min(Voc_array))/(max(Voc_array)-min(Voc_array));
% F = @(coef,xdata)coef(1)*exp(xdata*coef(2));
% monoexp_pendent_fit = fitnlm(Voc_array, monoexp_tau_array, 'Weigth', weigth);
% pendent = monoexp_pendent_fit.Coefficients.Estimate(2);

%% save results
TPV_struct.sol_name = symstructs{2, 1};
TPV_struct.Voc = Voc_array;
TPV_struct.monoexp_tau = monoexp_tau_array;
TPV_struct.monoexp_lin = monoexp_lin_array;
TPV_struct.biexp_tau_fast = biexp_tau_fast_array;
TPV_struct.biexp_tau_slow = biexp_tau_slow_array;
TPV_struct.biexp_lin_fast = biexp_lin_fast_array;
TPV_struct.biexp_lin_slow = biexp_lin_slow_array;

%% plot results
TPV_full_analysis(TPV_struct);
    
% make the figures appear, all at the end of the script
set(0, 'DefaultFigureVisible', 'on');
figHandles = findall(groot, 'Type', 'figure');
set(figHandles(:), 'visible', 'on')

%------------- END OF CODE --------------
