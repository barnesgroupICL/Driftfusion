function ISwave_struct = ISwave_full_exec(ssol_i_eq, ssol_i_light, startInt, endInt, Int_points, startFreq, endFreq, Freq_points, deltaV, BC, frozen_ions, frozen_freecharges, calcJi, parallelize, save_solutions, save_result)
% do Impedance Spectroscopy (IS) in a range of background light intensities
% applying an oscillating voltage

ssol_i_Int = ssol_i_light;
output_name = inputname(1);

% in case ions should be frozen during IS, set the mobility to zero after stabilization at new light intensity
original_p = ssol_i_Int.params;

% decrease annoiance by figures popping up
ssol_i_Int.params.figson = 0;
ssol_i_eq.params.figson = 0;
% don't display figures until the end of the script, as they steal the focus
% taken from https://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
set(0, 'DefaultFigureVisible', 'off');

tmesh_type = 1; % linear time mesh, 'cause here we're not interested in a decay, but in a oscillation going on during the whole simulated timespan

% number of complete oscillation periods to simulate
% the current looks reproducible already after few oscillations, this could be set in an automatic way
periods = 10;

% for having a meaningful output from verifyStabilization, here use a
% number of tpoints which is 1 + a multiple of 4 * periods
tpoints = 1 + 10 * 4 * periods; % gets redefined by changeLight, so re-setting is needed

% default pdepe tolerance is 1e-3, for having an accurate phase from the
% fitting, improving the tollerance is useful
RelTol = 1e-4;

% define light intensity values, always decreasing
Int_array = logspace(log10(startInt), log10(endInt), Int_points);

% in certain cases (e.g. no ions moving) dark takes too long to reach
% stability, so ssol_i_eq is used
Int_array = [Int_array, 0];
% anyway dark gives strange results

Int_array = unique(Int_array); % remove duplicates, useful just in case dark is selected as the only illumination condition from input options
Int_array = sort(Int_array, 'descend'); % from high intensity to low, needed by changelight in case also dark is included, otherwise logspace from zero to something breaks

% define frequency values
Freq_array = logspace(log10(startFreq), log10(endFreq), Freq_points);

% pre allocate arrays filling them with zeros
Voc_array = zeros(length(Int_array));
tmax_matrix = zeros(length(Int_array), length(Freq_array));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do a serie of IS measurements
disp('ISwave_full_exec - Doing the IS at various light intensities');
% set zero time for stabilization tmax, so changeLight will determine a good one
changeLight_tmax = 0;

for i = 1:length(Int_array)
    if frozen_ions
        ssol_i_Int.params.mui = original_p.mui; % reset default value of ion mobility in case frozen_ions option was used
    end
    if frozen_freecharges
        sol_i_Int.params.mue_i = original_p.mue_i; % electron mobility in intrinsic
        sol_i_Int.params.muh_i = original_p.muh_i; % hole mobility in intrinsic
        sol_i_Int.params.mue_p = original_p.mue_p;
        sol_i_Int.params.muh_p = original_p.muh_p;
        sol_i_Int.params.mue_n = original_p.mue_n;
        sol_i_Int.params.muh_n = original_p.muh_n;
    end
    if Int_array(i) % in dark use the ssol_i_eq solution, needed for no ions case where changeLight cannot reach dark
        ssol_i_Int = changeLight(ssol_i_Int, Int_array(i), changeLight_tmax); % decrease light intensity, on the symmetrical solution
        changeLight_tmax = ssol_i_Int.params.tmax / 2; % time to use for next iteration
    else
        ssol_i_Int = ssol_i_eq;
    end
    [sol_i_Int, Voc_array(i)] = asymmetricize(ssol_i_Int, BC); % normal BC 1 should work, also BC 2 can be employed
    if frozen_ions
        sol_i_Int.params.mui = 0; % if frozen_ions, freezing ions
    end
    if frozen_freecharges
        sol_i_Int.params.mue_i = 0; % electron mobility in intrinsic
        sol_i_Int.params.muh_i = 0; % hole mobility in intrinsic
        sol_i_Int.params.mue_p = 0;
        sol_i_Int.params.muh_p = 0;
        sol_i_Int.params.mue_n = 0;
        sol_i_Int.params.muh_n = 0;
        sol_i_Int.params.Int = 0;
    end 
    parfor (j = 1:length(Freq_array), parforArg)
        tempRelTol = RelTol; % convert RelTol variable to a temporary variable, as suggested for parallel loops
        sol_i_Int_ISwave = ISwave_single_exec(sol_i_Int, BC, Voc_array(i), deltaV, Freq_array(j), periods, tmesh_type, tpoints, calcJi, tempRelTol); % do IS
        if save_solutions && ~parallelize % assignin cannot be used in a parallel loop, so single solutions cannot be saved
            sol_i_Int_ISwave.params.figson = 1; % re-enable figures by default when using the saved solution, that were disabled above
            sol_name = matlab.lang.makeValidName([output_name '_Int_' num2str(Int_array(i)) '_Freq_' num2str(Freq_array(j)) '_ISwave']);
            assignin('base', sol_name, sol_i_Int_ISwave);
        end
        [fit_coeff, fit_idrift_coeff, ~, ~, ~, ~, ~, ~] = ISwave_single_analysis(sol_i_Int_ISwave, parallelize); % extract parameters and do plot
        % if phase is small or negative, double check increasing accuracy of the solver
        if fit_coeff(3) < 0.03
            disp('ISwave_full_exec - Fitted phase is very small or negative, double checking with higher solver accuracy')
            tempRelTol = tempRelTol / 100;
            sol_i_Int_ISwave = ISwave_single_exec(sol_i_Int, BC, Voc_array(i), deltaV, Freq_array(j), periods, tmesh_type, tpoints, calcJi, tempRelTol); % do IS
            [fit_coeff, fit_idrift_coeff, ~, ~, ~, ~, ~, ~] = ISwave_single_analysis(sol_i_Int_ISwave, parallelize); % repeat analysis on new solution
        end
        % if the phase is negative even with the new accuracy, check again
        if fit_coeff(3) < 0.003
            disp('ISwave_full_exec - Fitted phase is extremely small, increasing solver accuracy again')
            tempRelTol = tempRelTol / 100;
            sol_i_Int_ISwave = ISwave_single_exec(sol_i_Int, BC, Voc_array(i), deltaV, Freq_array(j), periods, tmesh_type, tpoints, calcJi, tempRelTol); % do IS
            [fit_coeff, fit_idrift_coeff, ~, ~, ~, ~, ~, ~] = ISwave_single_analysis(sol_i_Int_ISwave, parallelize); % repeat analysis on new solution
        end
        J_bias(i, j) = fit_coeff(1); % not really that useful
        J_amp(i, j) = fit_coeff(2);
        J_phase(i, j) = fit_coeff(3);
        J_i_bias(i, j) = fit_idrift_coeff(1);
        J_i_amp(i, j) = fit_idrift_coeff(2);
        J_i_phase(i, j) = fit_idrift_coeff(3);

        % as the number of periods is fixed, there's no need for tmax to be
        % a matrix, but this could change, so it's a matrix
        tmax_matrix(i,j) = sol_i_Int_ISwave.params.tmax;
    end
end

sun_index = find(Int_array == 1); % could used for plotting... maybe...

% even if here the frequency is always the same for each illumination, it
% is not the case for ISstep, and the solution has to be more similar in
% order to be used by the same IS_full_analysis_vsfrequency script
Freq_matrix = repmat(Freq_array, length(Int_array), 1);

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
pulsatance_matrix = 2 * pi * repmat(Freq_array, length(Int_array), 1);
% the capacitance is the imaginary part of 1/(pulsatance*complex_impedance)
% or can be obtained in the same way with Joutphase/(pulsatance*deltaV)
cap = sin(J_phase) ./ (pulsatance_matrix .* impedance_abs);

impedance_i_abs = -deltaV ./ J_i_amp; % J_amp is in amperes
impedance_i_re = impedance_i_abs .* cos(J_i_phase); % this is the resistance
impedance_i_im = impedance_i_abs .* sin(J_i_phase);
cap_idrift = sin(J_i_phase) ./ (pulsatance_matrix .* impedance_i_abs);

% save results, this struct is similar to ISstep_struct in terms of fields,
% but the columns and rows in the fields can be different
ISwave_struct.sol_name = inputname(1);
ISwave_struct.Voc = Voc_array;
ISwave_struct.periods = periods;
ISwave_struct.Freq = Freq_matrix;
ISwave_struct.tmesh_type = tmesh_type;
ISwave_struct.tpoints = tpoints;
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
    assignin('base', ['ISwave_' output_name], ISwave_struct);
end

% IS_full_analysis_vsvoltage(ISwave_struct); % there's no such script yet
IS_full_analysis_vsfrequency(ISwave_struct);
ISwave_full_analysis_nyquist(ISwave_struct);

% make the figures appear, all at the end of the script
set(0, 'DefaultFigureVisible', 'on');
figHandles = findall(groot, 'Type', 'figure');
set(figHandles(:), 'visible', 'on')
