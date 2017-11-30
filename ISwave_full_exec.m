function ISwave_struct = ISwave_full_exec(ssol_i_eq, ssol_i_light, startInt, endInt, Int_points, startFreq, endFreq, Freq_points, deltaV, BC, calcJi, save_solutions, save_results)
% do Impedance Spectroscopy (IS) in a range of background light intensities
% applying an oscillating voltage

ssol_i_Int = ssol_i_light;
% saved_mui = ssol_i_Int.params.mui;

% decrease annoiance by figures popping up
ssol_i_Int.params.figson = 0;
ssol_i_eq.params.figson = 0;
% don't display figures until the end of the script, as they steal the focus
% taken from https://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
set(0, 'DefaultFigureVisible', 'off');

periods = 4; % the current looks reproducible already after few oscillations, this could be set in an automatic way

% for having a meaningful output from verifyStabilization, here use a
% number of tpoints which is 1 + a multiple of 4 * periods
tpoints = 1 + 15 * 4 * periods; % gets redefined by changeLight, so re-setting is needed

% define light intensity values, always decreasing
Int_array = logspace(log10(startInt), log10(endInt), Int_points);

% in certain cases (e.g. no ions moving) dark takes too long to reach
% stability, so ssol_i_eq is used
Int_array = [Int_array, 0];

Int_array = [Int_array, 1]; % add 1 sun, for plotting the 1 sun Voc
Int_array = unique(Int_array); % remove duplicates
Int_array = sort(Int_array, 'descend'); % from high intensity to low

% define frequency values
Freq_array = logspace(log10(startFreq), log10(endFreq), Freq_points);

% pre allocate arrays filling them with zeros
Voc_array = zeros(length(Int_array));
tmax_matrix = zeros(length(Int_array), length(Freq_array));

% do a serie of IS measurements
disp('ISwave_full_exec - Doing the IS at various light intensities');
% set zero time for stabilization tmax, so changeLight will determine a good one
changeLight_tmax = 0;
for i = 1:length(Int_array)
    if Int_array(i) % in dark use the ssol_i_eq solution, needed for no ions case where changeLight cannot reach dark
        ssol_i_Int = changeLight(ssol_i_Int, Int_array(i), changeLight_tmax); % decrease light intensity, on the symmetrical solution
        changeLight_tmax = ssol_i_Int.params.tmax / 2; % time to use for next iteration
    else
        ssol_i_Int = ssol_i_eq;
    end
    [sol_i_Int, Voc_array(i)] = asymmetricize(ssol_i_Int, BC); % normal BC 1 should work, also BC 2 can be employed
    for j = 1:length(Freq_array)
        sol_i_Int_ISwave = ISwave_single_exec(sol_i_Int, BC, deltaV, Freq_array(j), periods, tpoints, calcJi); % do IS
        if save_solutions
            sol_name = matlab.lang.makeValidName([inputname(1) '_Int_' num2str(Int_array(i)) '_Freq_' num2str(Freq_array(j)) '_ISwave']);
            assignin('base', sol_name, sol_i_Int_ISwave);
        end
        % as the number of periods is fixed, there's no need for tmax to be
        % a matrix, but this could change, so it's a matrix
        tmax_matrix(i,j) = sol_i_Int_ISwave.params.tmax;
        [fit_coeff, ~, ~, ~, ~] = ISwave_single_analysis(sol_i_Int_ISwave); % extract parameters and do plot

    end
end

sun_index = find(Int_array == 1); % could used for plotting... maybe...





% save results, this struct is similar to ISstep_struct in terms of fields,
% but the columns and rows in the fields can be different
ISwave_struct.sol_name = inputname(1);
ISwave_struct.Voc = Voc_array;
ISwave_struct.sun_index = sun_index;
ISwave_struct.tmesh_type = sol_i_Int_ISwave.params.tmesh_type;
ISwave_struct.periods = periods;
ISwave_struct.tpoints = tpoints;
ISwave_struct.Int = Int_array;
ISwave_struct.BC = BC;
ISwave_struct.deltaV = deltaV;
ISwave_struct.J_bias = fit_coeff(1);
ISwave_struct.J_amp = fit_coeff(2);
ISwave_struct.J_phase = fit_coeff(3);
ISwave_struct.Freq = Freq_array;
ISwave_struct.cap = ;
ISwave_struct.resist = ;

if save_results
    assignin('base', ['ISwave_' inputname(1)], ISwave_struct);
end

IS_full_analysis_vsvoltage(ISwave_struct);
IS_full_analysis_vsfrequency(ISwave_struct);

% make the figures appear, all at the end of the script
set(0, 'DefaultFigureVisible', 'on');
figHandles = findall(groot, 'Type', 'figure');
set(figHandles(:), 'visible', 'on')
