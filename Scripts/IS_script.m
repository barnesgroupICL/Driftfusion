function IS_results = IS_script(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, demodulation, do_graphics)
%IS_SCRIPT - Simulates impedance spectroscopy at various frequencies on many provided solutions
% Alternatively, a cell containing various structures can be provided.
% For example, if a cell generated using genIntStructs is provided, the
% impedance at various light intensities can be compared.
% Or, if the provided cell has been generated using genVappStructs,
% solutions with different background voltage bias can be compared.
%
% Syntax:  IS_results = IS_script(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, demodulation, do_graphics)
%
% Inputs:
%   STRUCTS - can be a cell structure containing either symmetrical or
%     asymmetrical structures. Various background light intensities can be
%     provided using genIntStructs and various applied voltages
%     using genVappStructs.
%     Otherwise it can be a single struct as created by DF.
%   STARTFREQ - highest frequency limit
%   ENDFREQ - lowest frequency limit
%   FREQ_POINTS - number of points to simulate between STARTFREQ and
%     ENDFREQ
%   DELTAV - voltage oscillation amplitude in volts, one mV should be
%     enough but larger values can be employed for reducing the noise or
%     having a large perturbation simulation
%   FROZEN_IONS - logical, after stabilization sets the mobility of
%     ionic defects to zero to simulate the impedance with effectively
%     frozen ions
%   DEMODULATION - logical, determines which method to use for extracting
%     phase and amplitude of the current. If false, always uses fitting via
%     IS_EA_ana_fit, if true uses demodulation via
%     IS_EA_ana_demodulation. Anyway if the obtained phase is weird,
%     fit will be used automatically for confirming the result
%   DO_GRAPHICS - logical, whether to graph the individual solutions and
%     the overall graphics
%
% Outputs:
%   IS_RESULTS - a struct containing the most important results of the simulation
%
% Example:
%   IS_oc = IS_script(genIntStructsVoc(soleq.ion, 1, 1e-3, 7, true), 1e9, 1e-2, 56, 2e-3, false, true, true)
%     calculate starting from open circuit conditions on 8 different illumination intensities including dark, do not freeze ions, use a half peak to peak
%     voltage oscillation amplitude of 2 mV, on 23 points from frequencies of 1 GHz to
%     0.01 Hz and plot all graphics, including the ones for each solution
%   IS_oc_frozenions = IS_script(genIntStructsVoc(soleq.ion, 1, 1e-3, 7, true), 1e9, 1e-2, 56, 2e-3, true, true, true)
%     as above but freezing ions during voltage oscillation
%   IS_sc = IS_script(genIntStructs(soleq.ion, 1, 1e-3, 7, true), 1e9, 1e-2, 56, 2e-3, false, true, true)
%     as the first example but starting from short circuit conditions
%   IS_vapp = IS_script(genVappStructs(soleq.ion, 0:0.2:1, true), 1e9, 1e-2, 56, 2e-3, false, true, true)
%     as the first example but in dark and applying various voltages
%   IS_dark = IS_script(soleq.ion, 1e9, 1e-2, 56, 2e-3, false, true, true)
%     as above but using only dark solution with no bias voltage
%
% Other m-files required: doIS_EA, IS_ana_plot, IS_script_ana_nyquist, IS_script_ana_impedance, 
%   IS_script_ana_phase, dfana, stabilize
% Subfunctions: none
% MAT-files required: none
%
% See also genIntStructs, df, doIS_EA, IS_script_ana_nyquist, IS_ana_plot, IS_script_ana_impedance, IS_script_ana_phase.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

% in case a single struct is given in input, convert it to a cell structure
% with just one cell
if length(structs(:, 1)) == 1 % if the input is a single structure instead of a cell with structures
    structs_temp = cell(2, 1);
    structs_temp{1, 1} = structs;
    structs_temp{2, 1} = inputname(1);
    structs = structs_temp;
end

% don't display figures until the end of the script, as they steal the
% focus and are very annoying
% taken from https://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
% if do_graphics
%     set(0, 'DefaultFigureVisible', 'off');
% end

% number of complete oscillation periods to simulate
% the current looks reproducible already after few oscillations, this could be set in an automatic way
% this number should be above 20 for having good phase estimation in dark
% solutions via IS_EA_ana_demodulation
periods = 20;

% for having a meaningful output from verifyStabilization, here use a
% number of tpoints which is 1 + a multiple of 4 * periods
tpoints_per_period = 10 * 4; % gets redefined by changeLight, so re-setting is needed

% default pdepe tolerance is 1e-3, for having an accurate phase from the
% fitting, improving the tollerance is useful
RelTol = 1e-5;

% define frequency values
Freq_array = logspace(log10(startFreq), log10(endFreq), Freq_points);

%% pre allocate arrays filling them with zeros
Vdc_array = zeros(length(structs(1, :)), 1);
Int_array = Vdc_array;
tmax_matrix = zeros(length(structs(1, :)), length(Freq_array));
Jtot_bias = tmax_matrix;
Jtot_amp = tmax_matrix;
Jtot_phase = tmax_matrix;
ion_disp_bias = tmax_matrix;
ion_disp_amp = tmax_matrix;
ion_disp_phase = tmax_matrix;
cat_disp_bias = tmax_matrix;
cat_disp_amp = tmax_matrix;
cat_disp_phase = tmax_matrix;
ani_disp_bias = tmax_matrix;
ani_disp_amp = tmax_matrix;
ani_disp_phase = tmax_matrix;
r_bias = tmax_matrix;
r_amp = tmax_matrix;
r_phase = tmax_matrix;
np_dt_bias = tmax_matrix;
np_dt_amp = tmax_matrix;
np_dt_phase = tmax_matrix;

%% do a serie of IS measurements

disp([mfilename ' - Doing the IS']);
for i = 1:length(structs(1, :))
    s = structs{1, i};
    Int_array(i) = s.par.int1;
    % decrease annoiance by figures popping up
    s.par.figson = 0;
    % ensuring stability of starting solution
    s = stabilize(s);
    % calculate currently applied DC voltage, as defined in pinana
    Vapp_arr = dfana.calcVapp(s);
    Vdc_array(i) = Vapp_arr(end);
    Vdc_temp = Vdc_array(i); % convert from array to var for avoiding complaints from parfor
    % in case the simulation without moving ions is requested, freeze them
    if frozen_ions
        s.par.mobseti = false;
    end
    % if Parallel Computing Toolbox is not available, the following line
    % will work as a normal for cycle
    parfor (j = 1:length(Freq_array), Inf)
        tempRelTol = RelTol; % convert RelTol variable to a temporary variable, as suggested for parallel loops
        s_IS = doIS_EA(s, deltaV,...
            Freq_array(j), periods, tpoints_per_period, 0.5, tempRelTol); % do IS
        % set IS_single_analysis minimal_mode to true as under
        % parallelization graphics for single solutions cannot be created
        coeff = IS_ana_plot(s_IS, false, demodulation);
        % a phase close to 90 degrees can be indicated as it was -90 degree
        % by demodulation in case the RelTol was not enough

        % if phase is small or negative, double check increasing accuracy of the solver
        if coeff.Jtot(3) < 0.006 || coeff.Jtot(3) > pi/2 - 0.006
            disp([mfilename ' - Int: ' num2str(s.par.int1) '; Vdc: ' num2str(Vdc_temp) ' V; Freq: ' num2str(Freq_array(j)) ' Hz; Phase is ' num2str(rad2deg(coeff.Jtot(3))) ' degrees, increasing solver accuracy and calculating again'])
            % decrease tollerance
            tempRelTol = tempRelTol / 100;
            % start from the oscillating solution last point, better starting point
            s_IS = doIS_EA(s_IS,...
                deltaV, Freq_array(j), periods, tpoints_per_period, 0.5, tempRelTol); % do IS
            % repeat analysis on new solution
            coeff = IS_ana_plot(s_IS, false, demodulation);
        end
        % if phase is negative or bigger than pi/2, it could be a failure of demodulation or a real thing,
        % for confirming that is a real thing use the alternative fitting method without repeating the simulation
        if coeff.Jtot(3) < 0 || coeff.Jtot(3) > pi/2
            disp([mfilename ' - Int: ' num2str(s.par.int1) '; Vdc: ' num2str(Vdc_temp) ' V; Freq: ' num2str(Freq_array(j)) ' Hz; Phase is weird: ' num2str(rad2deg(coeff.Jtot(3))) ' degrees, confirming using alternative fitting method'])
            % in case demodulation was being used, use fitting instead, just in case the weird result was due to
            % demodulation failure
            coeff = IS_ana_plot(s_IS, false, ~demodulation);
            disp([mfilename ' - Int: ' num2str(s.par.int1) '; Vdc: ' num2str(Vdc_temp) ' V; Freq: ' num2str(Freq_array(j)) ' Hz; Phase from alternative fitting method is: ' num2str(rad2deg(coeff.Jtot(3))) ' degrees'])
        end
        % if phase is still negative or more than pi/2, check again increasing accuracy
        if coeff.Jtot(3) < 0 || abs(coeff.Jtot(3)) > pi/2
            disp([mfilename ' - Int: ' num2str(s.par.int1) '; Vdc: ' num2str(Vdc_temp) ' V; Freq: ' num2str(Freq_array(j)) ' Hz; Phase is still weird: ' num2str(rad2deg(coeff.Jtot(3))) ' degrees, increasing solver accuracy and calculating again'])
            tempRelTol = tempRelTol / 100;
            % start from the oscillating solution, better starting point
            s_IS = doIS_EA(s_IS,...
                deltaV, Freq_array(j), periods, tpoints_per_period, 0.5, tempRelTol); % do IS
            coeff = IS_ana_plot(s_IS, false, demodulation);
        end
        % save values
        Jtot_bias(i, j) = coeff.Jtot(1);
        Jtot_amp(i, j) = coeff.Jtot(2);
        Jtot_phase(i, j) = coeff.Jtot(3);
        ion_disp_bias(i, j) = coeff.ion_disp(1);
        ion_disp_amp(i, j) = coeff.ion_disp(2);
        ion_disp_phase(i, j) = coeff.ion_disp(3);
        cat_disp_bias(i, j) = coeff.cat_disp(1);
        cat_disp_amp(i, j) = coeff.cat_disp(2);
        cat_disp_phase(i, j) = coeff.cat_disp(3);
        if isfield(coeff, 'ani_disp')
            ani_disp_bias(i, j) = coeff.ani_disp(1);
            ani_disp_amp(i, j) = coeff.ani_disp(2);
            ani_disp_phase(i, j) = coeff.ani_disp(3);
        end
        r_bias(i, j) = coeff.r(1);
        r_amp(i, j) = coeff.r(2);
        r_phase(i, j) = coeff.r(3);
        np_dt_bias(i, j) = coeff.np_dt(1);
        np_dt_amp(i, j) = coeff.np_dt(2);
        np_dt_phase(i, j) = coeff.np_dt(3);

        
        % as the number of periods is fixed, there's no need for tmax to be
        % a matrix, but this could change in future code...
        tmax_matrix(i,j) = s_IS.par.tmax;
    end
end

%% calculate apparent capacity and impedance

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
impedance_abs = deltaV ./ Jtot_amp; % J_amp is in amperes
% the components of the impedance gets calculated with the phase from the
% current-voltage "delay"
% impedance phase is minus current phase, so -J_phase
impedance_re = impedance_abs .* cos(-Jtot_phase); % this is the resistance
impedance_im = impedance_abs .* sin(-Jtot_phase);
pulsatance_matrix = 2 * pi * repmat(Freq_array, length(structs(1, :)), 1);
% the capacitance is the imaginary part of 1/(pulsatance*complex_impedance)
% or can be obtained in the same way with Joutphase/(pulsatance*deltaV)
cap = sin(Jtot_phase) ./ (pulsatance_matrix .* impedance_abs);

%% impedance due to ionic displacement current

impedance_ion_disp_abs = deltaV ./ ion_disp_amp; % J_amp is in amperes
% impedance phase is minus current phase, so -J_i_phase
impedance_ion_disp_re = impedance_ion_disp_abs .* cos(-ion_disp_phase); % this is the resistance
impedance_ion_disp_im = impedance_ion_disp_abs .* sin(-ion_disp_phase);
cap_ion_disp = sin(ion_disp_phase) ./ (pulsatance_matrix .* impedance_ion_disp_abs);

impedance_cat_disp_abs = deltaV ./ cat_disp_amp; % J_amp is in amperes
% impedance phase is minus current phase, so -J_i_phase
impedance_cat_disp_re = impedance_cat_disp_abs .* cos(-cat_disp_phase); % this is the resistance
impedance_cat_disp_im = impedance_ion_disp_abs .* sin(-cat_disp_phase);
cap_cat_disp = sin(cat_disp_phase) ./ (pulsatance_matrix .* impedance_cat_disp_abs);

impedance_ani_disp_abs = deltaV ./ ani_disp_amp; % J_amp is in amperes
% impedance phase is minus current phase, so -J_i_phase
impedance_ani_disp_re = impedance_ani_disp_abs .* cos(-ani_disp_phase); % this is the resistance
impedance_ani_disp_im = impedance_ani_disp_abs .* sin(-ani_disp_phase);
cap_ani_disp = sin(ani_disp_phase) ./ (pulsatance_matrix .* impedance_ani_disp_abs);

%% impedance due to recombination current

impedance_r_abs = deltaV ./ r_amp; % J_amp is in amperes
% impedance phase is minus current phase, so -J_U_phase
impedance_r_re = impedance_r_abs .* cos(-r_phase); % this is the resistance
impedance_r_im = impedance_r_abs .* sin(-r_phase);
cap_r = sin(r_phase) ./ (pulsatance_matrix .* impedance_r_abs);

%% impedance due to accumulating current

impedance_np_dt_abs = deltaV ./ np_dt_amp; % J_amp is in amperes
impedance_np_dt_re = impedance_np_dt_abs .* cos(-np_dt_phase); % this is the resistance
impedance_np_dt_im = impedance_np_dt_abs .* sin(-np_dt_phase);
cap_np_dt = sin(np_dt_phase) ./ (pulsatance_matrix .* impedance_np_dt_abs);

%% save results

IS_results.sol_name = structs{2, 1};
IS_results.Vdc = Vdc_array;
IS_results.periods = periods;
IS_results.Freq = Freq_matrix;
IS_results.tpoints = 1 + tpoints_per_period * periods;
IS_results.tmax = tmax_matrix;
IS_results.Int = Int_array;
IS_results.deltaV = deltaV;
IS_results.sun_index = sun_index;
IS_results.Jtot_bias = Jtot_bias;
IS_results.Jtot_amp = Jtot_amp;
IS_results.Jtot_phase = Jtot_phase;
IS_results.ion_disp_bias = ion_disp_bias;
IS_results.ion_disp_amp = ion_disp_amp;
IS_results.ion_disp_phase = ion_disp_phase;
IS_results.cat_disp_bias = cat_disp_bias;
IS_results.cat_disp_amp = cat_disp_amp;
IS_results.cat_disp_phase = cat_disp_phase;
IS_results.ani_disp_bias = ani_disp_bias;
IS_results.ani_disp_amp = ani_disp_amp;
IS_results.ani_disp_phase = ani_disp_phase;
IS_results.r_bias = r_bias;
IS_results.r_amp = r_amp;
IS_results.r_phase = r_phase;
IS_results.np_dt_bias = np_dt_bias;
IS_results.np_dt_amp = np_dt_amp;
IS_results.np_dt_phase = np_dt_phase;
IS_results.cap = cap;
IS_results.impedance_abs = impedance_abs;
IS_results.impedance_im = impedance_im;
IS_results.impedance_re = impedance_re;
IS_results.cap_ion_disp = cap_ion_disp;
IS_results.impedance_ion_disp_abs = impedance_ion_disp_abs;
IS_results.impedance_ion_disp_im = impedance_ion_disp_im;
IS_results.impedance_ion_disp_re = impedance_ion_disp_re;
IS_results.cap_cat_disp = cap_cat_disp;
IS_results.impedance_cat_disp_abs = impedance_cat_disp_abs;
IS_results.impedance_cat_disp_im = impedance_cat_disp_im;
IS_results.impedance_cat_disp_re = impedance_cat_disp_re;
IS_results.cap_ani_disp = cap_ani_disp;
IS_results.impedance_ani_disp_abs = impedance_ani_disp_abs;
IS_results.impedance_ani_disp_im = impedance_ani_disp_im;
IS_results.impedance_ani_disp_re = impedance_ani_disp_re;
IS_results.cap_r = cap_r;
IS_results.impedance_r_abs = impedance_r_abs;
IS_results.impedance_r_im = impedance_r_im;
IS_results.impedance_r_re = impedance_r_re;
IS_results.cap_np_dt = cap_np_dt;
IS_results.impedance_np_dt_abs = impedance_np_dt_abs;
IS_results.impedance_np_dt_im = impedance_np_dt_im;
IS_results.impedance_np_dt_re = impedance_np_dt_re;

%% plot results

if do_graphics
    IS_script_plot_phase(IS_results);
    IS_script_plot_impedance(IS_results);
    IS_script_plot_nyquist(IS_results);
end

% make the figures appear, all at the end of the script
% set(0, 'DefaultFigureVisible', 'on');
% figHandles = findall(groot, 'Type', 'figure');
% set(figHandles(:), 'visible', 'on')

%------------- END OF CODE --------------
