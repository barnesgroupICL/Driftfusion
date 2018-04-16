function ISwave_results = ISwave_full_exec(structs, startFreq, endFreq, Freq_points, deltaV, BC, frozen_ions, do_graphics)
%ISWAVE_FULL_EXEC - Do Impedance Spectroscopy approximated applying an
% oscillating voltage (ISwave) in a range of background light intensities
%
% Syntax:  ISwave_results = ISwave_full_exec(structs, startFreq, endFreq, Freq_points, deltaV, BC, frozen_ions, do_graphics)
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
%   BC - boundary conditions indicating if the contacts are selective, see
%     PINDRIFT
%   FROZEN_IONS - logical, after stabilization sets the mobility of
%     ionic defects to zero
%   DO_GRAPHICS - logical, whether to graph the individual solutions and
%     the overall graphics
%
% Outputs:
%   ISWAVE_RESULTS - a struct containing the most important results of the simulation
%
% Example:
%   ISwave_full_exec(genIntStructs(ssol_i_eq, ssol_i_light, 100, 1e-7, 4), 1e9, 1e-2, 23, 1e-3, 1, true, false, true)
%     calculate also with dark background, do not freeze ions, use a
%     voltage oscillation amplitude of 1 mV, on 23 points from frequencies of 1 GHz to
%     0.01 Hz, with selective contacts, without calculating ionic current,
%     without parallelization
%   ISwave_full_exec(genIntStructs(ssol_i_eq, ssol_i_light, 100, 1e-7, 4), 1e9, 1e-2, 23, 1e-3, 1, true, false, true)
%     as above but freezing ions during voltage oscillation
%   ISwave_full_exec(genIntStructs(ssol_i_eq, ssol_i_light, 100, 1e-7, 4), 1e9, 1e-2, 23, 1e-3, 1, true, true, true)
%     calculate ionic current in the middle of the intrinsic
%   ISwave_full_exec(ssol_i_light_BC2, 1e9, 1e-2, 23, 1e-3, 2, true, false, false)
%     use non perfectly selective contacts (BC = 2)
%
% Other m-files required: asymmetricize, ISwave_EA_single_exec,
%   ISwave_single_analysis, ISwave_full_analysis_nyquist,
%   IS_full_analysis_impedance, ISwave_EA_full_analysis_phase, pinAna
% Subfunctions: none
% MAT-files required: none
%
% See also genIntStructs, pindrift, ISwave_EA_single_exec, ISwave_full_analysis_nyquist, ISwave_single_analysis, IS_full_analysis_impedance, ISwave_EA_full_analysis_phase.

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
    structs_temp = cell(2, 1);
    structs_temp{1, 1} = structs;
    structs_temp{2, 1} = inputname(1);
    structs = structs_temp;
end

% don't display figures until the end of the script, as they steal the focus
% taken from https://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
if do_graphics
    set(0, 'DefaultFigureVisible', 'off');
end

% which method to use for extracting phase and amplitude of the current
% if false, always uses fitting, if true uses demodulation multiplying the current
% by sin waves. Anyway if the obtained phase is werid, fit will be used
demodulation = true;

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

disp([mfilename ' - Doing the IS at various light intensities']);
for i = 1:length(structs(1, :))
    struct = structs{1, i};
    Int_array(i) = struct.p.Int;
    % decrease annoiance by figures popping up
    struct.p.figson = 0;
    if struct.p.OC % in case the solution is symmetric, break it in halves
        asymstruct_Int = asymmetricize(struct, BC); % normal BC 1 should work, also BC 2 can be employed
    else
        asymstruct_Int = struct;
    end
    [~, ~, ~, Efn, Efp, ~] = pinAna(asymstruct_Int);
    Vdc_array(i) = Efn(end, end) - Efp(end, 1);
    if frozen_ions
        asymstruct_Int.p.mui = 0; % if frozen_ions option is set, freezing ions
    end
    % if Parallel Computing Toolbox is not available, the following line
    % will work as a normal for cycle
    parfor (j = 1:length(Freq_array), Inf)
        tempRelTol = RelTol; % convert RelTol variable to a temporary variable, as suggested for parallel loops
        asymstruct_ISwave = ISwave_EA_single_exec(asymstruct_Int, BC, deltaV,...
            Freq_array(j), periods, tpoints_per_period, true, false, tempRelTol); % do IS
        % set ISwave_single_analysis minimal_mode to true as under
        % parallelization graphics for single solutions cannot be created
        [n_coeff, i_coeff, U_coeff, dQ_coeff] = ISwave_single_analysis(asymstruct_ISwave, true, demodulation);
        % if phase is small or negative, double check increasing accuracy of the solver
        % a phase close to 90 degrees can be indicated as it was -90 degree
        % by the demodulation, the fitting way does not have this problem
        if n_coeff(3) < 0.006 || n_coeff(3) > pi/2 - 0.006
            disp([mfilename ' - Freq: ' num2str(Freq_array(j)) '; Fitted phase is ' num2str(rad2deg(n_coeff(3))) ' degrees, it is extremely small or close to pi/2 or out of 0-pi/2 range, increasing solver accuracy and calculate again'])
            tempRelTol = tempRelTol / 100;
            % start from the oscillating solution, better starting point
            asymstruct_ISwave = ISwave_EA_single_exec(asymstruct_ISwave, BC,...
                deltaV, Freq_array(j), periods, tpoints_per_period, true, false, tempRelTol); % do IS
            % set ISwave_single_analysis minimal_mode is true if parallelize is true
            % repeat analysis on new solution
            [n_coeff, i_coeff, U_coeff, dQ_coeff] = ISwave_single_analysis(asymstruct_ISwave, true, demodulation);
        end
        % if phase is still negative or bigger than pi/2, likely is demodulation that is
        % failing (no idea why), use safer fitting method without repeating
        % the simulation
        if n_coeff(3) < 0 || n_coeff(3) > pi/2
            disp([mfilename ' - Freq: ' num2str(Freq_array(j)) '; Phase from demodulation is weird: ' num2str(rad2deg(n_coeff(3))) ' degrees, confirming using fitting'])
            % use fitting
            [n_coeff, i_coeff, U_coeff, dQ_coeff] = ISwave_single_analysis(asymstruct_ISwave, true, false);
            disp([mfilename ' - Freq: ' num2str(Freq_array(j)) '; Phase from fitting is: ' num2str(rad2deg(n_coeff(3))) ' degrees'])
        end
        % if phase is still negative or more than pi/2, check again increasing accuracy
        if n_coeff(3) < 0 || abs(n_coeff(3)) > pi/2
            disp([mfilename ' - Freq: ' num2str(Freq_array(j)) '; Fitted phase is ' num2str(rad2deg(n_coeff(3))) ' degrees, it is out of 0-pi/2 range, increasing solver accuracy and calculate again'])
            tempRelTol = tempRelTol / 100;
            % start from the oscillating solution, better starting point
            asymstruct_ISwave = ISwave_EA_single_exec(asymstruct_ISwave, BC,...
                deltaV, Freq_array(j), periods, tpoints_per_period, true, false, tempRelTol); % do IS
            % set ISwave_single_analysis minimal_mode is true if parallelize is true
            % repeat analysis on new solution
            % use fitting
            [n_coeff, i_coeff, U_coeff, dQ_coeff] = ISwave_single_analysis(asymstruct_ISwave, true, false);
        end
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
        % a matrix, but this could change, so it's a matrix
        tmax_matrix(i,j) = asymstruct_ISwave.p.tmax;
    end
end

%% calculate apparent capacity 

sun_index = find(Int_array == 1); % could used for plotting... maybe...

% even if here the frequency is always the same for each illumination, it
% is not the case for ISstep, and the solution has to be more similar in
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

impedance_i_abs = deltaV ./ J_i_amp; % J_amp is in amperes
% impedance phase is minus current phase, so -J_i_phase
impedance_i_re = impedance_i_abs .* cos(-J_i_phase); % this is the resistance
impedance_i_im = impedance_i_abs .* sin(-J_i_phase);
cap_idrift = sin(J_i_phase) ./ (pulsatance_matrix .* impedance_i_abs);

impedance_U_abs = deltaV ./ J_U_amp; % J_amp is in amperes
% impedance phase is minus current phase, so -J_U_phase
impedance_U_re = impedance_U_abs .* cos(-J_U_phase); % this is the resistance
impedance_U_im = impedance_U_abs .* sin(-J_U_phase);
cap_U = sin(J_U_phase) ./ (pulsatance_matrix .* impedance_U_abs);

impedance_dQ_abs = deltaV ./ dQ_amp; % J_amp is in amperes
impedance_dQ_re = impedance_dQ_abs .* cos(-dQ_phase); % this is the resistance
impedance_dQ_im = impedance_dQ_abs .* sin(-dQ_phase);
cap_dQ = sin(dQ_phase) ./ (pulsatance_matrix .* impedance_dQ_abs);

%% save results

% this struct is similar to ISstep_struct in terms of fields,
% but the columns and rows in the fields can be different
ISwave_results.sol_name = structs{2, 1};
ISwave_results.Vdc = Vdc_array;
ISwave_results.periods = periods;
ISwave_results.Freq = Freq_matrix;
ISwave_results.tpoints = 1 + tpoints_per_period * periods;
ISwave_results.tmax = tmax_matrix;
ISwave_results.Int = Int_array;
ISwave_results.BC = BC;
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
    ISwave_EA_full_analysis_phase(ISwave_results);
    IS_full_analysis_impedance(ISwave_results);
    ISwave_full_analysis_nyquist(ISwave_results);
end

% make the figures appear, all at the end of the script
set(0, 'DefaultFigureVisible', 'on');
figHandles = findall(groot, 'Type', 'figure');
set(figHandles(:), 'visible', 'on')

%------------- END OF CODE --------------
