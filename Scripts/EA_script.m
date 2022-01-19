function EA_results = EA_script(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, do_graphics)
%EA_SCRIPT - Do ElectroAbsorption (Stark spectroscopy, observing electric field in the absorber applying an oscillating voltage) in a range of background light intensities
% this file is heavily based on the equivalent script for Impedance
% Spectroscopy: IS_script
% Similarly to the impedance simulation performed with
% IS_script, this routine applies an oscillating voltage
% on all the provided solutions.
% One relevant difference is the fact that the current calculation is
% disabled, making the analysis faster.
%
% Syntax:  EA_results = EA_script(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, do_graphics)
%
% Inputs:
%   STRUCTS - can be a cell structure containing structs at various background
%     light intensities. This can be generated using genIntStructs.
%     Otherwise it can be a single struct as created by DF.
%   STARTFREQ - higher frequency limit
%   ENDFREQ - lower frequency limit
%   FREQ_POINTS - number of points to simulate between STARTFREQ and
%     ENDFREQ
%   DELTAV - voltage oscillation amplitude in volts, one mV should be enough
%   FROZEN_IONS - logical, after stabilization sets the mobility of
%     ionic defects to zero
%   DO_GRAPHICS - logical, whether to graph the individual solutions and
%     the overall graphics
%
% Outputs:
%   EA_RESULTS - a struct containing the most important results of the simulation
%
% Example:
%   EA_ill = EA_script(genIntStructs(soleq.ion, 100, 1e-7, 4, true), 1e9, 1e-2, 23, 1e-3, false, true)
%     calculate at illuminations from 100 suns to dark, do not freeze ions, use a
%     voltage oscillation amplitude of 1 mV, on 23 points from frequencies of 1 GHz to
%     0.01 Hz, plotting the results
%   EA_ill_frozenions = EA_script(genIntStructs(soleq.ion, 100, 1e-7, 4, true), 1e9, 1e-2, 23, 1e-3, true, true)
%     as above but freezing ions during voltage oscillation
%   EA_dark = EA_script(soleq.ion, 1e9, 1e-2, 23, 1e-3, true, true)
%     as above but only with no background illumination
%
% Other m-files required: doIS_EA, EA_ana_plot, EA_script_ana_phase, dfana, EA_script_ana_Efield
% Subfunctions: none
% MAT-files required: none
%
% See also genIntStructs, df, doIS_EA, EA_script_ana_phase, EA_ana_plot, EA_script_ana_Efield.

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

% don't display figures until the end of the script, as they steal the focus
% taken from https://stackoverflow.com/questions/8488758/inhibit-matlab-window-focus-stealing
% if do_graphics
%     set(0, 'DefaultFigureVisible', 'off');
% end

% which method to use for extracting phase and amplitude of the current.
% If false: uses fitting. If true: uses demodulation.
demodulation = true;

% number of complete oscillation periods to simulate
% the current looks reproducible already after few oscillations, this could be set in an automatic way
% this number should be above 20 for having good phase estimation in dark
% solutions via IS_EA_ana_demodulation
periods = 20;

% for having a meaningful output from verifyStabilization, here use a
% number of tpoints which is 1 + a multiple of 4 * periods
tpoints_per_period = 20 * 4; % gets redefined by changeLight, so re-setting is needed

% default pdepe tolerance is 1e-3, for having an accurate phase from the
% fitting, improving the tollerance is useful
RelTol = 1e-4;

% define frequency values
Freq_array = logspace(log10(startFreq), log10(endFreq), Freq_points);

%% pre allocate arrays filling them with zeros
Vdc_array = zeros(length(structs(1, :)), 1);
Int_array = Vdc_array;
tmax_matrix = zeros(length(structs(1, :)), length(Freq_array));

%% do a serie of EA simulations

disp([mfilename ' - Doing the EA at various light intensities']);
for i = 1:length(structs(1, :))
    s = structs{1, i};
    Int_array(i) = s.par.int1;
    % decrease annoiance by figures popping up
    s.par.figson = 0;
    % ensuring stability of starting solution
    s = stabilize(s);
    Vapp_arr = dfana.calcVapp(s);
    Vdc_array(i) = Vapp_arr(end);
    if frozen_ions
        s.par.mobseti = 0; % if frozen_ions option is set, freezing ions
    end
    % if Parallel Computing Toolbox is not available, the following line
    % will work as a normal for cycle
    parfor (j = 1:length(Freq_array), Inf)
        tempRelTol = RelTol; % convert RelTol variable to a temporary variable, as suggested for parallel loops
        struct_EA = doIS_EA(s, deltaV,...
            Freq_array(j), periods, tpoints_per_period, 0.5, tempRelTol); % do EA
        % set EA_single_analysis do_graphics to false as under
        % parallelization graphics for single solutions cannot be created
        % set local_field to false for not calculating the field for each
        % position in the intrinsic
        coeff = EA_ana_plot(struct_EA, false, false, demodulation, missing);
        AC_ExDC_E_bias(i, j) =     coeff.h1_mean(1);
        AC_ExDC_E_amp(i, j) =      coeff.h1_mean(2);
        AC_ExDC_E_phase(i, j) =    coeff.h1_mean(3);
        AC_ExDC_E_i_bias(i, j) =   coeff.i_h1_mean(1);
        AC_ExDC_E_i_amp(i, j) =    coeff.i_h1_mean(2);
        AC_ExDC_E_i_phase(i, j) =  coeff.i_h1_mean(3);
        AC_Efield2_bias(i, j) =    coeff.h2_mean(1);
        AC_Efield2_amp(i, j) =     coeff.h2_mean(2);
        AC_Efield2_phase(i, j) =   coeff.h2_mean(3);
        AC_Efield2_i_bias(i, j) =  coeff.i_h2_mean(1);
        AC_Efield2_i_amp(i, j) =   coeff.i_h2_mean(2);
        AC_Efield2_i_phase(i, j) = coeff.i_h2_mean(3);
        AC_Efield_amp_squared_mean(i, j) = coeff.ac_amp_squared_mean;
        
        % as the number of periods is fixed, there's no need for tmax to be
        % a matrix, but this could change, so it's a matrix
        tmax_matrix(i,j) = struct_EA.par.tmax;
    end
end

%% calculate

sun_index = find(Int_array == 1); % could used for plotting... maybe...

% even if here the frequency is always the same for each illumination, it
% is not the case for ISstep, and the solution has to be more similar in
% order to be used by the same scripts
Freq_matrix = repmat(Freq_array, length(structs(1, :)), 1);

%% save results

% this struct is similar to ISstep_struct and IS_struct in terms of fields,
% but the columns and rows in the fields can be different
EA_results.sol_name = structs{2, 1};
EA_results.Vdc = Vdc_array;
EA_results.periods = periods;
EA_results.Freq = Freq_matrix;
EA_results.tpoints = 1 + tpoints_per_period * periods;
EA_results.tmax = tmax_matrix;
EA_results.Int = Int_array;
EA_results.deltaV = deltaV;
EA_results.sun_index = sun_index;
EA_results.AC_ExDC_E_bias = AC_ExDC_E_bias;
EA_results.AC_ExDC_E_amp = AC_ExDC_E_amp;
EA_results.AC_ExDC_E_phase = AC_ExDC_E_phase;
EA_results.AC_ExDC_E_i_bias = AC_ExDC_E_i_bias;
EA_results.AC_ExDC_E_i_amp = AC_ExDC_E_i_amp;
EA_results.AC_ExDC_E_i_phase = AC_ExDC_E_i_phase;
EA_results.AC_Efield2_bias = AC_Efield2_bias;
EA_results.AC_Efield2_amp = AC_Efield2_amp;
EA_results.AC_Efield2_phase = AC_Efield2_phase;
EA_results.AC_Efield2_i_bias = AC_Efield2_i_bias;
EA_results.AC_Efield2_i_amp = AC_Efield2_i_amp;
EA_results.AC_Efield2_i_phase = AC_Efield2_i_phase;
EA_results.AC_Efield_amp_squared_mean = AC_Efield_amp_squared_mean;

%% plot results

if do_graphics
    EA_script_plot_phase(EA_results);
    EA_script_plot_Efield(EA_results);
end

% make the figures appear, all at the end of the script
% set(0, 'DefaultFigureVisible', 'on');
% figHandles = findall(groot, 'Type', 'figure');
% set(figHandles(:), 'visible', 'on')

%------------- END OF CODE --------------
