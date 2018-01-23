function ISwave_struct = ISwave_full_exec(symstructs, startFreq, endFreq, Freq_points, deltaV, BC, frozen_ions, calcJi, parallelize, save_solutions, save_result)
% do Impedance Spectroscopy (IS) in a range of background light intensities
% applying an oscillating voltage

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

% number of complete oscillation periods to simulate
% the current looks reproducible already after few oscillations, this could be set in an automatic way
periods = 10;

% for having a meaningful output from verifyStabilization, here use a
% number of tpoints which is 1 + a multiple of 4 * periods
tpoints_per_period = 10 * 4; % gets redefined by changeLight, so re-setting is needed

% default pdepe tolerance is 1e-3, for having an accurate phase from the
% fitting, improving the tollerance is useful
RelTol = 1e-4;

% define frequency values
Freq_array = logspace(log10(startFreq), log10(endFreq), Freq_points);

% pre allocate arrays filling them with zeros
Voc_array = zeros(length(symstructs(1, :)));
Int_array = Voc_array;
tmax_matrix = zeros(length(symstructs(1, :)), length(Freq_array));
J_bias = tmax_matrix;
J_amp = tmax_matrix;
J_phase = tmax_matrix;
J_i_bias = tmax_matrix;
J_i_amp = tmax_matrix;
J_i_phase = tmax_matrix;

% switch on or off the parallelization calculation of solution, when true
% less pictures gets created and the solutions does not get saved in main
% workspace
if parallelize
  parforArg = Inf;
else
  parforArg = 0;
end

%% do a serie of IS measurements

disp([mfilename ' - Doing the IS at various light intensities']);
for i = 1:length(symstructs(1, :))
    Int_array(i) = symstructs{1, i}.params.Int;
    % decrease annoiance by figures popping up
    symstructs{1, i}.params.figson = 0;
    [asymstruct_Int, Voc_array(i)] = asymmetricize(symstructs{1, i}, BC); % normal BC 1 should work, also BC 2 can be employed
    if frozen_ions
        asymstruct_Int.params.mui = 0; % if frozen_ions option is set, freezing ions
    end
    parfor (j = 1:length(Freq_array), parforArg)
        tempRelTol = RelTol; % convert RelTol variable to a temporary variable, as suggested for parallel loops
        asymstruct_ISwave = ISwave_single_exec(asymstruct_Int, BC, deltaV, Freq_array(j), periods, tpoints_per_period, calcJi, tempRelTol); % do IS
        if save_solutions && ~parallelize % assignin cannot be used in a parallel loop, so single solutions cannot be saved
            sol_name = matlab.lang.makeValidName([symstructs{2, i} '_Freq_' num2str(Freq_array(j)) '_ISwave']);
            asymstruct_ISwave.params.figson = 1; % re-enable figures by default when using the saved solution, that were disabled above
            assignin('base', sol_name, asymstruct_ISwave);
        end
        % set ISwave_single_analysis minimal_mode is true if parallelize is true
        [fit_coeff, fit_idrift_coeff, ~, ~, ~, ~, ~, ~] = ISwave_single_analysis(asymstruct_ISwave, parallelize); % extract parameters and do plot
        % if phase is small or negative, double check increasing accuracy of the solver
        if fit_coeff(3) < 0.03
            disp([mfilename ' - Fitted phase is very small or negative, double checking with higher solver accuracy'])
            tempRelTol = tempRelTol / 100;
            asymstruct_ISwave = ISwave_single_exec(asymstruct_Int, BC, deltaV, Freq_array(j), periods, tpoints_per_period, calcJi, tempRelTol); % do IS
            % set ISwave_single_analysis minimal_mode is true if parallelize is true
            [fit_coeff, fit_idrift_coeff, ~, ~, ~, ~, ~, ~] = ISwave_single_analysis(asymstruct_ISwave, parallelize); % repeat analysis on new solution
        end
        % if the phase is negative even with the new accuracy, check again
        if fit_coeff(3) < 0.003
            disp([mfilename ' - Fitted phase is extremely small, increasing solver accuracy again'])
            tempRelTol = tempRelTol / 100;
            asymstruct_ISwave = ISwave_single_exec(asymstruct_Int, BC, deltaV, Freq_array(j), periods, tpoints_per_period, calcJi, tempRelTol); % do IS
            % set ISwave_single_analysis minimal_mode is true if parallelize is true
            [fit_coeff, fit_idrift_coeff, ~, ~, ~, ~, ~, ~] = ISwave_single_analysis(asymstruct_ISwave, parallelize); % repeat analysis on new solution
        end
        J_bias(i, j) = fit_coeff(1); % not really that useful
        J_amp(i, j) = fit_coeff(2);
        J_phase(i, j) = fit_coeff(3);
        J_i_bias(i, j) = fit_idrift_coeff(1);
        J_i_amp(i, j) = fit_idrift_coeff(2);
        J_i_phase(i, j) = fit_idrift_coeff(3);

        % as the number of periods is fixed, there's no need for tmax to be
        % a matrix, but this could change, so it's a matrix
        tmax_matrix(i,j) = asymstruct_ISwave.params.tmax;
    end
end

sun_index = find(Int_array == 1); % could used for plotting... maybe...

% even if here the frequency is always the same for each illumination, it
% is not the case for ISstep, and the solution has to be more similar in
% order to be used by the same IS_full_analysis_vsfrequency script
Freq_matrix = repmat(Freq_array, length(symstructs(1, :)), 1);

% deltaV is a scalar, J_amp and J_phase are matrices
% as the current of MPP is defined as positive in the model, we expect that
% with a positive deltaV we have a negative J_amp (J_amp is forced to be negative actually)

% the absolute value of impedance has to be taken from the absolute values
% of oscillation of voltage and of current
impedance_abs = -deltaV ./ J_amp; % J_amp is in amperes
% the components of the impedance gets calculated with the phase from the
% current-voltage "delay"
impedance_re = impedance_abs .* cos(J_phase); % this is the resistance
impedance_im = impedance_abs .* sin(J_phase);
pulsatance_matrix = 2 * pi * repmat(Freq_array, length(symstructs(1, :)), 1);
% the capacitance is the imaginary part of 1/(pulsatance*complex_impedance)
% or can be obtained in the same way with Joutphase/(pulsatance*deltaV)
cap = sin(J_phase) ./ (pulsatance_matrix .* impedance_abs);

impedance_i_abs = -deltaV ./ J_i_amp; % J_amp is in amperes
impedance_i_re = impedance_i_abs .* cos(J_i_phase); % this is the resistance
impedance_i_im = impedance_i_abs .* sin(J_i_phase);
cap_idrift = sin(J_i_phase) ./ (pulsatance_matrix .* impedance_i_abs);

% save results, this struct is similar to ISstep_struct in terms of fields,
% but the columns and rows in the fields can be different
ISwave_struct.sol_name = symstructs{2, 1};
ISwave_struct.Voc = Voc_array;
ISwave_struct.periods = periods;
ISwave_struct.Freq = Freq_matrix;
ISwave_struct.tpoints = 1 + tpoints_per_period * periods;
ISwave_struct.tmax = tmax_matrix;
ISwave_struct.Int = Int_array;
ISwave_struct.BC = BC;
ISwave_struct.deltaV = deltaV;
ISwave_struct.sun_index = sun_index;
ISwave_struct.J_bias = J_bias;
ISwave_struct.J_amp = J_amp;
ISwave_struct.J_phase = J_phase;
ISwave_struct.J_i_bias = J_i_bias;
ISwave_struct.J_i_amp = J_i_amp;
ISwave_struct.J_i_phase = J_i_phase;
ISwave_struct.cap = cap;
ISwave_struct.impedance_abs = impedance_abs;
ISwave_struct.impedance_im = impedance_im;
ISwave_struct.impedance_re = impedance_re;
ISwave_struct.cap_idrift = cap_idrift;
ISwave_struct.impedance_i_abs = impedance_i_abs;
ISwave_struct.impedance_i_im = impedance_i_im;
ISwave_struct.impedance_i_re = impedance_i_re;
if save_result
    assignin('base', ['ISwave_' symstructs{2, 1}], ISwave_struct);
end

IS_full_analysis_vsfrequency(ISwave_struct);
ISwave_full_analysis_nyquist(ISwave_struct);

% make the figures appear, all at the end of the script
set(0, 'DefaultFigureVisible', 'on');
figHandles = findall(groot, 'Type', 'figure');
set(figHandles(:), 'visible', 'on')
