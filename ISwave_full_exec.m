function ISwave_struct = ISwave_full_exec(ssol_i_eq, ssol_i_light, startInt, endInt, Int_points, startFreq, endFreq, Freq_points, deltaV, BC, frozen_ions, calcJi, parallelize, save_solutions, save_results)
% do Impedance Spectroscopy (IS) in a range of background light intensities
% applying an oscillating voltage

ssol_i_Int = ssol_i_light;
output_name = inputname(1);

% in case ions should be frozen during IS, set the mobility to zero after stabilization at new light intensity
saved_mui = ssol_i_Int.params.mui;
if frozen_ions
    new_mui = 0;
else
    new_mui = saved_mui;
end

% decrease annoiance by figures popping up
ssol_i_Int.params.figson = 0;
ssol_i_eq.params.figson = 0;
% don't display figures until the end of the script, as they steal the focus
% taken from https://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
set(0, 'DefaultFigureVisible', 'off');

tmesh_type = 1;
periods = 10; % the current looks reproducible already after few oscillations, this could be set in an automatic way

% for having a meaningful output from verifyStabilization, here use a
% number of tpoints which is 1 + a multiple of 4 * periods
tpoints = 1 + 10 * 4 * periods; % gets redefined by changeLight, so re-setting is needed

% default pdepe tolerance is 1e-3, in case of small phase, BetterRelTol is
% used instead
BetterRelTol = 1e-5;

% define light intensity values, always decreasing
Int_array = logspace(log10(startInt), log10(endInt), Int_points);

% in certain cases (e.g. no ions moving) dark takes too long to reach
% stability, so ssol_i_eq is used
Int_array = [Int_array, 0];
% anyway dark gives strange results

% Int_array = unique(Int_array); % remove duplicates
% Int_array = sort(Int_array, 'descend'); % from high intensity to low

% define frequency values
Freq_array = logspace(log10(startFreq), log10(endFreq), Freq_points);

% pre allocate arrays filling them with zeros
Voc_array = zeros(length(Int_array));
tmax_matrix = zeros(length(Int_array), length(Freq_array));
J_bias = tmax_matrix;
J_amp = tmax_matrix;
J_phase = tmax_matrix;
J_idrift_bias = tmax_matrix;
J_idrift_amp = tmax_matrix;
J_idrift_phase = tmax_matrix;

% switch on or off the parallelization calculation of solution, when true
% less pictures gets created and the solutions does not get saved in main
% workspace
if parallelize
  parforArg = Inf;
else
  parforArg = 0;
end

% do a serie of IS measurements
disp('ISwave_full_exec - Doing the IS at various light intensities');
% set zero time for stabilization tmax, so changeLight will determine a good one
changeLight_tmax = 0;

tic
for i = 1:length(Int_array)
%     RelTol = 1e-3; % reset default value
    ssol_i_Int.params.mui = saved_mui; % reset default value of ion mobility in case frozen_ions option was used
    if Int_array(i) % in dark use the ssol_i_eq solution, needed for no ions case where changeLight cannot reach dark
        ssol_i_Int = changeLight(ssol_i_Int, Int_array(i), changeLight_tmax); % decrease light intensity, on the symmetrical solution
        changeLight_tmax = ssol_i_Int.params.tmax / 2; % time to use for next iteration
    else
        ssol_i_Int = ssol_i_eq;
    end
    [sol_i_Int, Voc_array(i)] = asymmetricize(ssol_i_Int, BC); % normal BC 1 should work, also BC 2 can be employed
    sol_i_Int.params.mui = new_mui; % if frozen_ions, freezing ions
    parfor (j = 1:length(Freq_array), parforArg)
        sol_i_Int_ISwave = ISwave_single_exec(sol_i_Int, BC, Voc_array(i), deltaV, Freq_array(j), periods, tmesh_type, tpoints, calcJi, BetterRelTol); % do IS
        [fit_coeff, fit_idrift_coeff, ~, ~, ~, ~, ~, ~] = ISwave_single_analysis(sol_i_Int_ISwave); % extract parameters and do plot
        % if phase is small or negative, double check increasing accuracy of the solver
%         if RelTol ~= BetterRelTol && fit_coeff(3) < 0.05
%             disp('ISwave_full_exec - Fitted phase is very small or negative, double checking with higher solver accuracy')
%             RelTol = BetterRelTol; % this stays set until a new light intensity is tested
%             % here I start from the oscillating solution sol_i_Int_ISwave, 
%             % which is more risky, beware
%             sol_i_Int_ISwave = ISwave_single_exec(sol_i_Int, BC, Voc_array(i), deltaV, Freq_array(j), periods, tpoints, calcJi, RelTol); % do IS
%             [fit_coeff, ~, ~, ~, ~] = ISwave_single_analysis(sol_i_Int_ISwave); % repeat analysis on new solution
%         end
%         if RelTol == BetterRelTol && fit_coeff(3) > 0.05
%             disp('ISwave_full_exec - Fitted phase is not small anymore, decreasing solver accuracy')
%             RelTol = 1e-3; % reset default value if the phase is not small anymore
%         end
        J_bias(i, j) = fit_coeff(1); % not really that useful
        J_amp(i, j) = fit_coeff(2);
        J_phase(i, j) = fit_coeff(3);
        J_idrift_bias(i, j) = fit_idrift_coeff(1);
        J_idrift_amp(i, j) = fit_idrift_coeff(2);
        J_idrift_phase(i, j) = fit_idrift_coeff(3);
        if save_solutions && ~parallelize
            sol_name = matlab.lang.makeValidName([output_name '_Int_' num2str(Int_array(i)) '_Freq_' num2str(Freq_array(j)) '_ISwave']);
            assignin('base', sol_name, sol_i_Int_ISwave);
        end
        % as the number of periods is fixed, there's no need for tmax to be
        % a matrix, but this could change, so it's a matrix
        tmax_matrix(i,j) = sol_i_Int_ISwave.params.tmax;
    end
end
toc

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
impedance = -deltaV ./ J_amp; % J_amp is in amperes
% the components of the impedance gets calculated with the phase from the
% current-voltage "delay"
impedance_re = impedance .* cos(J_phase); % this is the resistance
impedance_im = impedance .* sin(J_phase);
pulsatance_matrix = 2 * pi * repmat(Freq_array, length(Int_array), 1);
% the capacitance is the imaginary part of 1/(pulsatance*complex_impedance)
% or can be obtained in the same way with Joutphase/(pulsatance*deltaV)
cap = sin(J_phase) ./ (pulsatance_matrix .* impedance);

impedance_idrift = -deltaV ./ J_idrift_amp; % J_amp is in amperes
impedance_idrift_re = impedance_idrift .* cos(J_idrift_phase); % this is the resistance
impedance_idrift_im = impedance_idrift .* sin(J_idrift_phase);
cap_idrift = sin(J_idrift_phase) ./ (pulsatance_matrix .* impedance_idrift);

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
ISwave_struct.J_idrift_bias = J_idrift_bias;
ISwave_struct.J_idrift_amp = J_idrift_amp;
ISwave_struct.J_idrift_phase = J_idrift_phase;
ISwave_struct.cap = cap;
ISwave_struct.impedance_im = impedance_im;
ISwave_struct.impedance_re = impedance_re;
ISwave_struct.cap_idrift = cap_idrift;
ISwave_struct.impedance_idrift_im = impedance_idrift_im;
ISwave_struct.impedance_idrift_re = impedance_idrift_re;
if save_results
    assignin('base', ['ISwave_' output_name], ISwave_struct);
end

% IS_full_analysis_vsvoltage(ISwave_struct); % there's no such script yet
IS_full_analysis_vsfrequency(ISwave_struct);
ISwave_full_analysis_nyquist(ISwave_struct);

% make the figures appear, all at the end of the script
set(0, 'DefaultFigureVisible', 'on');
figHandles = findall(groot, 'Type', 'figure');
set(figHandles(:), 'visible', 'on')
