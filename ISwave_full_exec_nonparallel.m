function ISwave_results = ISwave_full_exec_nonparallel(structs, startFreq, endFreq, Freq_points, deltaV, sequential, frozen_ions, demodulation, do_graphics, save_solutions)
%ISWAVE_FULL_EXEC_NONPARALLEL - Do Impedance Spectroscopy approximated applying an oscillating voltage (ISwave) in a range of background light intensities.
% Getting rid of annoying parfor
%
% Syntax:  ISwave_results = ISwave_full_exec_nonparallel(structs, startFreq, endFreq, Freq_points, deltaV, sequential, frozen_ions, demodulation, do_graphics, save_solutions)
%
% Inputs:
%   STRUCTS - can be a cell structure containing structs at various background
%     light intensities. This can be generated using genIntStructs.
%     Otherwise it can be a single struct as created by PINDRIFT.
%   STARTFREQ - higher frequency limit
%   ENDFREQ - lower frequency limit
%   FREQ_POINTS - number of points to simulate between STARTFREQ and
%     ENDFREQ
%   DELTAV - voltage oscillation amplitude in volts, one mV should be enough
%   SEQUENTIAL - logical, if true do not check that the oscillating solution reached a
%     (oscillating) stabilization, instead take always the first solution
%     and use its final timepoint as the starting point of the next
%     frequency simulation. This can be useful when it's known that the
%     starting solution is not stabilized and a realistic simulation of the
%     solution evolving during the measurement is wanted.
%   FROZEN_IONS - logical, after stabilization sets the mobility of
%     ionic defects to zero
%   DEMODULATION - which method to use for extracting phase and amplitude of the current
%     if false, always uses fitting, if true uses demodulation multiplying the current
%     by sin waves. Anyway if the obtained phase is werid, fit will be used
%     automatically for confirming the result
%   DO_GRAPHICS - logical, whether to graph the individual solutions and
%     the overall graphics
%   SAVE_SOLUTIONS - is a logic defining if to assing in volatile base
%     workspace the calulated solutions of single ISstep perturbations
%
% Outputs:
%   ISWAVE_RESULTS - a struct containing the most important results of the simulation
%
% Example:
%   ISwave_oc = ISwave_full_exec_nonparallel(genIntStructs(ssol_i_eq_SR, 1, 1e-3, 7, true), 1e9, 1e-2, 56, 2e-3, false, false, true, true, true)
%     calculate on 8 different illumination intensities including dark,
%     on 23 points from frequencies of 1 GHz to 0.01 Hz,
%     use a half peak to peak voltage oscillation amplitude of 2 mV,
%     start each simulation from a stabilized solution,
%     do not freeze ions,
%     using demodulation instead of fitting,
%     plot all graphics, including the ones for each solution,
%     save to base workspace the single solutions
%   ISwave_oc_sequential = ISwave_full_exec_nonparallel(genIntStructs(ssol_i_eq_SR, 1, 1e-3, 7, true), 1e9, 1e-2, 56, 2e-3, true, false, true, true, true)
%     as above but starting each simulation from the last point of the
%     previous simulation, without reaching a stable solution
%   ISwave_sc = ISwave_full_exec_nonparallel(genIntStructs(sol_i_eq_SR, 1, 1e-3, 7, true), 1e9, 1e-2, 56, 2e-3, false, false, true, true, true)
%     as the first example but starting from short circuit conditions
%   ISwave_vapp = ISwave_full_exec_nonparallel(genVappStructs(sol_i_eq_SR, 0:0.2:1), 1e9, 1e-2, 56, 2e-3, false, false, true, true, true)
%     as the first example but in dark and applying various voltages
%
% Other m-files required: asymmetricize, ISwave_EA_single_exec,
%   ISwave_single_analysis, ISwave_full_analysis_nyquist,
%   IS_full_analysis_impedance, ISwave_full_analysis_phase, pinana,
%   stabilize
% Subfunctions: none
% MAT-files required: none
%
% See also genIntStructs, pindrift, ISwave_EA_single_exec, ISwave_full_analysis_nyquist, ISwave_single_analysis, IS_full_analysis_impedance, ISwave_full_analysis_phase.

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

% don't display figures until the end of the script, as they steal the
% focus and being very annoying
% taken from https://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
if do_graphics
    set(0, 'DefaultFigureVisible', 'off');
end

% number of complete oscillation periods to simulate
% the current looks reproducible already after few oscillations, this could be set in an automatic way
% this number should be above 20 for having good phase estimation in dark
% solutions via ISwave_EA_single_demodulation
periods = 20;

% for having a meaningful output from verifyStabilization, here use a
% number of tpoints which is 1 + a multiple of 4 * periods
tpoints_per_period = 10 * 4; % gets redefined by changeLight, so re-setting is needed

% default pdepe tolerance is 1e-3, for having an accurate phase from the
% fitting, improving the tollerance is useful
RelTol = 1e-6;

% define frequency values
Freq_array = logspace(log10(startFreq), log10(endFreq), Freq_points);

%% pre allocate arrays filling them with zeros
Vdc_array = zeros(length(structs(1, :)), 1);
Int_array = Vdc_array;
tmax_matrix = zeros(length(structs(1, :)), length(Freq_array));
J_bias = tmax_matrix;
J_amp = tmax_matrix;
J_phase = tmax_matrix;
J_i_bias = tmax_matrix;
J_i_amp = tmax_matrix;
J_i_phase = tmax_matrix;
J_U_bias = tmax_matrix;
J_U_amp = tmax_matrix;
J_U_phase = tmax_matrix;
dQ_bias = tmax_matrix;
dQ_amp = tmax_matrix;
dQ_phase = tmax_matrix;

%% do a serie of IS measurements

disp([mfilename ' - Doing the IS']);
for i = 1:length(structs(1, :))
    struct = structs{1, i};
    Int_array(i) = struct.p.Int;
    % decrease annoiance by figures popping up
    struct.p.figson = 0;
    if struct.p.OC % in case the solution is symmetric, break it in halves
        asymstruct_Int = asymmetricize(struct);
    else % otherwise use it as it comes
        asymstruct_Int = struct;
    end
    % the asymmetricized solution could require some stabilization after
    % breaking (mostly when both high illumination and mobile ions are present)
    if ~sequential
        asymstruct_Int = stabilize(asymstruct_Int);
    end
    % calculate currently applied DC voltage, as defined in pinana
    [Vapp_arr, ~, ~] = pinana(asymstruct_Int);
    Vdc_array(i) = Vapp_arr(end);
    % in case the simulation without moving ions is requested, freeze them
    if frozen_ions
        asymstruct_Int.p.mui = 0;
    end
    % simulate first frequency with the stabilized solution
    asymstruct_start = asymstruct_Int;
    for j = 1:length(Freq_array)
        tempRelTol = RelTol; % reset RelTol variable
        asymstruct_ISwave = ISwave_EA_single_exec(asymstruct_start, deltaV,...
            Freq_array(j), periods, tpoints_per_period, ~sequential, false, tempRelTol); % do IS
        % set ISwave_single_analysis minimal_mode to false if graphics are
        % requested
        % extract parameters and do plot
        [n_coeff, i_coeff, U_coeff, dQ_coeff, ~] = ISwave_single_analysis(asymstruct_ISwave, ~do_graphics, demodulation);
        % a phase close to 90 degrees can be indicated as it was -90 degree
        % by demodulation in case RelTol was not enough

        % if phase is small or negative, double check increasing accuracy of the solver
        if n_coeff(3) < 0.006 || n_coeff(3) > pi/2 - 0.006
            disp([mfilename ' - Int: ' num2str(Int_array(i)) '; Vdc: ' num2str(Vdc_array(i)) ' V; Freq: ' num2str(Freq_array(j)) ' Hz; Phase is ' num2str(rad2deg(n_coeff(3))) ' degrees, increasing solver accuracy and calculating again'])
            % decrease tollerance
            tempRelTol = tempRelTol / 100;
            % if just the initial solution, non-stabilized, is requested, do
            % not start from oscillating solution
            if sequential
                asymstruct_temp = asymstruct_start; % strictly use the last point from previous cycle
            else
                asymstruct_temp = asymstruct_ISwave; % the last point of the oscillating solution, better starting point
            end
            asymstruct_ISwave = ISwave_EA_single_exec(asymstruct_temp,...
                deltaV, Freq_array(j), periods, tpoints_per_period, ~sequential, false, tempRelTol); % do IS
            % repeat analysis on new solution
            [n_coeff, i_coeff, U_coeff, dQ_coeff, ~] = ISwave_single_analysis(asymstruct_ISwave, ~do_graphics, demodulation);
        end
        % if phase is negative or bigger than pi/2, it could be a failure of demodulation or a real thing,
        % for confirming that is a real thing use the alternative fitting method without repeating the simulation
        if n_coeff(3) < 0 || n_coeff(3) > pi/2
            disp([mfilename ' - Int: ' num2str(Int_array(i)) '; Vdc: ' num2str(Vdc_array(i)) ' V; Freq: ' num2str(Freq_array(j)) ' Hz; Phase is weird: ' num2str(rad2deg(n_coeff(3))) ' degrees, confirming using alternative fitting method'])
            % in case demodulation was being used, use fitting instead, just in case the weird result was due to
            % demodulation failure

            [n_coeff, i_coeff, U_coeff, dQ_coeff, ~] = ISwave_single_analysis(asymstruct_ISwave, ~do_graphics, ~demodulation);
            disp([mfilename ' - Int: ' num2str(Int_array(i)) '; Vdc: ' num2str(Vdc_array(i)) ' V; Freq: ' num2str(Freq_array(j)) ' Hz; Phase from alternative fitting method is: ' num2str(rad2deg(n_coeff(3))) ' degrees'])
        end
        % if phase is still negative or more than pi/2, check again increasing accuracy
        if n_coeff(3) < 0 || abs(n_coeff(3)) > pi/2
            disp([mfilename ' - Int: ' num2str(Int_array(i)) '; Vdc: ' num2str(Vdc_array(i)) ' V; Freq: ' num2str(Freq_array(j)) ' Hz; Phase is still weird: ' num2str(rad2deg(n_coeff(3))) ' degrees, increasing solver accuracy and calculating again'])
            tempRelTol = tempRelTol / 100;
            % if just the initial solution, non-stabilized, is requested, do
            % not start from oscillating solution
            if sequential
                asymstruct_temp = asymstruct_start; % strictly use the last point from previous cycle
            else
                asymstruct_temp = asymstruct_ISwave; % the oscillating solution, better starting point
            end
            asymstruct_ISwave = ISwave_EA_single_exec(asymstruct_temp,...
                deltaV, Freq_array(j), periods, tpoints_per_period, ~sequential, false, tempRelTol); % do IS
            % repeat analysis on new solution
            [n_coeff, i_coeff, U_coeff, dQ_coeff, ~] = ISwave_single_analysis(asymstruct_ISwave, ~do_graphics, demodulation);
        end
        % save values
        J_bias(i, j) = n_coeff(1);
        J_amp(i, j) = n_coeff(2);
        J_phase(i, j) = n_coeff(3);
        J_i_bias(i, j) = i_coeff(1);
        J_i_amp(i, j) = i_coeff(2);
        J_i_phase(i, j) = i_coeff(3);
        J_U_bias(i, j) = U_coeff(1);
        J_U_amp(i, j) = U_coeff(2);
        J_U_phase(i, j) = U_coeff(3);
        dQ_bias(i, j) = dQ_coeff(1);
        dQ_amp(i, j) = dQ_coeff(2);
        dQ_phase(i, j) = dQ_coeff(3);
        
        % as the number of periods is fixed, there's no need for tmax to be
        % a matrix, but this could change in future code...
        tmax_matrix(i,j) = asymstruct_ISwave.p.tmax;
        
        if save_solutions % assignin cannot be used in a parallel loop, so single solutions cannot be saved
            sol_name = matlab.lang.makeValidName([structs{2, i} '_Freq_' num2str(Freq_array(j)) '_ISwave']);
            asymstruct_ISwave.p.figson = 1; % re-enable figures by default when using the saved solution, that were disabled above
            assignin('base', sol_name, asymstruct_ISwave);
        end
        % in case the next simulation is run starting from the last time
        % point of this simulation, redefine the starting point
        % otherwise the originally provided solution is used again
        % in case of big perturbations, to use the last point of previous
        % solution would make sense even when not using sequential. This is
        % because the background DC profiles are quite different (because
        % of big perturbations) from the stable provided input
        if sequential
            asymstruct_start = asymstruct_ISwave;
        end
    end
end

%% calculate apparent capacity  and impedance

% save in an easy to access variable which solution refers to 1 sun
sun_index = find(Int_array == 1);

% even if here the frequency is always the same for each illumination, it
% is not the case for ISstep (still unpublished), and the solution has to be more similar in
% order to be used by the same scripts
Freq_matrix = repmat(Freq_array, length(structs(1, :)), 1);

% deltaV is a scalar, J_amp and J_phase are matrices
% as the current of MPP is defined as positive in the model, we expect that
% with a positive deltaV we have a negative J_amp (J_amp is forced to be negative actually)

% the absolute value of impedance has to be taken from the absolute values
% of oscillation of voltage and of current
impedance_abs = deltaV ./ J_amp; % J_amp is in amperes
% the components of the impedance gets calculated with the phase from the
% current-voltage "delay"
% impedance phase is minus current phase, so -J_phase
impedance_re = impedance_abs .* cos(-J_phase); % this is the resistance
impedance_im = impedance_abs .* sin(-J_phase);
pulsatance_matrix = 2 * pi * repmat(Freq_array, length(structs(1, :)), 1);
% the capacitance is the imaginary part of 1/(pulsatance*complex_impedance)
% or can be obtained in the same way with Joutphase/(pulsatance*deltaV)
cap = sin(J_phase) ./ (pulsatance_matrix .* impedance_abs);

%% impedance due to ionic displacement current

impedance_i_abs = deltaV ./ J_i_amp; % J_amp is in amperes
% impedance phase is minus current phase, so -J_i_phase
impedance_i_re = impedance_i_abs .* cos(-J_i_phase); % this is the resistance
impedance_i_im = impedance_i_abs .* sin(-J_i_phase);
cap_idrift = sin(J_i_phase) ./ (pulsatance_matrix .* impedance_i_abs);

%% impedance due to recombination current

impedance_U_abs = deltaV ./ J_U_amp; % J_amp is in amperes
% impedance phase is minus current phase, so -J_U_phase
impedance_U_re = impedance_U_abs .* cos(-J_U_phase); % this is the resistance
impedance_U_im = impedance_U_abs .* sin(-J_U_phase);
cap_U = sin(J_U_phase) ./ (pulsatance_matrix .* impedance_U_abs);

%% impedance due to accumulating current

impedance_dQ_abs = deltaV ./ dQ_amp; % J_amp is in amperes
impedance_dQ_re = impedance_dQ_abs .* cos(-dQ_phase); % this is the resistance
impedance_dQ_im = impedance_dQ_abs .* sin(-dQ_phase);
cap_dQ = sin(dQ_phase) ./ (pulsatance_matrix .* impedance_dQ_abs);

%% save results

ISwave_results.sol_name = structs{2, 1};
ISwave_results.Vdc = Vdc_array;
ISwave_results.periods = periods;
ISwave_results.Freq = Freq_matrix;
ISwave_results.tpoints = 1 + tpoints_per_period * periods;
ISwave_results.tmax = tmax_matrix;
ISwave_results.Int = Int_array;
ISwave_results.deltaV = deltaV;
ISwave_results.sun_index = sun_index;
ISwave_results.J_bias = J_bias;
ISwave_results.J_amp = J_amp;
ISwave_results.J_phase = J_phase;
ISwave_results.J_i_bias = J_i_bias;
ISwave_results.J_i_amp = J_i_amp;
ISwave_results.J_i_phase = J_i_phase;
ISwave_results.J_U_bias = J_U_bias;
ISwave_results.J_U_amp = J_U_amp;
ISwave_results.J_U_phase = J_U_phase;
ISwave_results.dQ_bias = dQ_bias;
ISwave_results.dQ_amp = dQ_amp;
ISwave_results.dQ_phase = dQ_phase;
ISwave_results.cap = cap;
ISwave_results.impedance_abs = impedance_abs;
ISwave_results.impedance_im = impedance_im;
ISwave_results.impedance_re = impedance_re;
ISwave_results.cap_idrift = cap_idrift;
ISwave_results.impedance_i_abs = impedance_i_abs;
ISwave_results.impedance_i_im = impedance_i_im;
ISwave_results.impedance_i_re = impedance_i_re;
ISwave_results.cap_U = cap_U;
ISwave_results.impedance_U_abs = impedance_U_abs;
ISwave_results.impedance_U_im = impedance_U_im;
ISwave_results.impedance_U_re = impedance_U_re;
ISwave_results.cap_dQ = cap_dQ;
ISwave_results.impedance_dQ_abs = impedance_dQ_abs;
ISwave_results.impedance_dQ_im = impedance_dQ_im;
ISwave_results.impedance_dQ_re = impedance_dQ_re;

%% plot results

if do_graphics
    ISwave_full_analysis_phase(ISwave_results);
    IS_full_analysis_impedance(ISwave_results);
    ISwave_full_analysis_nyquist(ISwave_results);
end

% make the figures appear, all at the end of the script
set(0, 'DefaultFigureVisible', 'on');
figHandles = findall(groot, 'Type', 'figure');
set(figHandles(:), 'visible', 'on')

%------------- END OF CODE --------------
