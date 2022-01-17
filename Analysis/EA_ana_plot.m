function coeff = EA_ana_plot(struct_IS_EA, do_graphics, local_field, demodulation, savefig_dir)
%EA_ANA_PLOT - Calculate ElectroAbsorbance (ElectroChromism, Stark spectroscopy) and phase using the same output as the impedance spectroscopy with oscillating voltage
% this file is heavily based on the equivalent script for Impedance
% Spectroscopy: IS_ana_plot
% The first and second harmonics are locally calculated and averaged over the device thickness.
%
% Syntax:  EA_ana_plot(struct_IS_EA, do_graphics, local_field, demodulation, savefig_dir)
%
% Inputs:
%   STRUCT_IS_EA - a struct where a solution has been perturbed by an
%     oscillating voltage, as generated from doIS_EA
%   DO_GRAPHICS - logical, when true graphics does not get created, useful when
%     launched under parallelization
%   LOCAL_FIELD - logical, whether to fit the position resolved data
%   DEMODULATION - logical, get phase via demodulation instead of using a fitting
%   SAVEFIG_DIR - string, set the directory where to save .fig and .png
%     image files for second harmonic versus position.
%     To save in the current directory set the dot directory ".".
%     To disable saving, set "missing" value (without quotes).
%
% Outputs:
%   coeff.h1_mean - an array of bias, amplitude and phase of the first
%     harmonic
%   coeff.h1 - a matrix of bias, amplitude and phase of the first
%     harmonic, each column is a space point
%   coeff.i_h1_mean - an array of bias, amplitude and phase of the first
%     harmonic for ionic profile
%   coeff.i_h1 - a matrix of bias, amplitude and phase of the first
%     harmonic for ionic profile, each column is a space point
%   coeff.h2_mean - an array of bias, amplitude and phase of the second
%     harmonic
%   coeff.h2 - a matrix of bias, amplitude and phase of the second
%     harmonic, each column is a space point
%   coeff.i_h2_mean - an array of bias, amplitude and phase of the second
%     harmonic for ionic profile
%   coeff.i_h2 - a matrix of bias, amplitude and phase of the second
%     harmonic for ionic profile, each column is a space point
%   coeff.ac_amp_squared_mean - the value of the average squared amplitude of the
%     AC electric field. It would be identical to the amplitude of the
%     second harmonic if the contributions from each spatial point was in
%     phase
%   coeff.ac - a matrix of bias, amplitude and phase of the AC electric
%     field, each column is a space point
%
% Example:
%   EA_ana_plot(doIS_EA(soleq.ion, 2e-3, 1e6, 20, 40, 0.5, 1e-4), true, true, true, missing)
%     do simulation and plot, do not save images to files
%   EA_ana_plot(sol_i_eq_SR_ea_800mV_10mHz, true, true, true, ".")
%     do plot, save images to current directory
%
% Other m-files required: IS_EA_ana_fit, IS_EA_ana_demodulation
% Subfunctions: none
% MAT-files required: none
%
% See also EA_script, doIS_EA, IS_EA_ana_demodulation, IS_EA_ana_fit.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

% shortcut
s = struct_IS_EA;

% Intrinsic points logical array
absorber_points = logical(s.par.dev.g0);
absorber_points_num = sum(absorber_points);

%% preallocate arrays and fallback values

% pre allocate and set result to zero in case next code block is
% skipped due to local_field option
coeff.h1_mean = NaN(3, 1);
coeff.h1 = NaN(3, absorber_points_num);
coeff.i_h1_mean = NaN(3, 1);
coeff.i_h1 = NaN(3, absorber_points_num);
coeff.h2_mean = NaN(3, 1);
coeff.h2 = NaN(3, absorber_points_num);
coeff.i_h2_mean = NaN(3, 1);
coeff.i_h2 = NaN(3, absorber_points_num);
coeff.ac_amp_squared_mean = NaN;
coeff.ac = NaN(3, absorber_points_num);

%% check the input and prepare environment

% verify if the simulation broke, in that case return NaNs
if size(s.u, 1) < s.par.tpoints
    return;
end

% round should not be needed here
% s.par.V_fun_arg(3) should be frequency
periods = round(s.par.tmax * s.par.V_fun_arg(3));

% here is critical that exactely an entire number of periods is provided to
% IS_EA_ana_demodulation
fit_t_index = round((s.par.tpoints - 1) * round(periods / 2) / periods) + 1;
fit_t = s.t(fit_t_index:end)';

% obtain params for Vapp profile for second armonic demodulation, the
% result has the same frequency and phase as squared Vapp
%Vapp2_params = s.par.V_fun_arg;
%Vapp2_params(3) = Vapp2_params(3) * 2; % double the frequency
%Vapp2_params(4) = Vapp2_params(4) - pi/2; % subtract 90 degrees

%% calculate E total (DC + AC)

% for obtaining an averaged E, I could just take the potential difference
% between the last point in the intrinsic and the first, but I need to
% calculate the E field for then squaring it anyway

% obtain averaged E field in the intrinsic, which is the absorber
tot_Efield_everywhere = dfana.calcF(s);
tot_Efield = tot_Efield_everywhere(:, absorber_points);
%gradient(s.u(:, absorber_points, 1), s.x(absorber_points), s.t); % I don't know why here s.t is needed as it's not being used

%% calculate E DC

% obtain averaged E field at the last time point in the intrinsic, which is the absorber
%DC_Efield = gradient(sDC.u(end, itype_points, 1), sDC.x(itype_points));
% the problem with the previous way is that in the middle of the material
% can easily result zero, and also the correct value to take is the average
% of the total field (AC + DC) which could (see the case for ions) be
% different from the DC value
DC_Efield = trapz(s.t(fit_t_index:end), tot_Efield(fit_t_index:end, :), 1) / (s.t(end) - s.t(fit_t_index));

%% calculate E AC

% tot_Efield is a matrix, DC_Efield is a vector
AC_Efield = tot_Efield - DC_Efield;
fit_AC_Efield = AC_Efield(fit_t_index:end, :);

%% calculate first harmonic, AC E field * DC E field

% the first harmonic is proportional to DC electric field multiplied by AC
% electric field

% AC_Efield is a matrix, DC_Efield is an array
AC_ExDC_E = AC_Efield .* DC_Efield;
fit_AC_ExDC_E = AC_ExDC_E(fit_t_index:end, :);
% average over absorber thickness, keep just time resolution
absorber_x = s.x(absorber_points);
absorber_thickness = absorber_x(end)-absorber_x(1);
AC_ExDC_E_mean = trapz(absorber_x, AC_ExDC_E, 2) / absorber_thickness;
fit_AC_ExDC_E_mean = AC_ExDC_E_mean(fit_t_index:end);

if demodulation
    coeff.h1_mean = IS_EA_ana_demodulation(fit_t, fit_AC_ExDC_E_mean, s.par.V_fun_type, s.par.V_fun_arg(3));
else
    coeff.h1_mean = IS_EA_ana_fit(fit_t, fit_AC_ExDC_E_mean, s.par.V_fun_type, s.par.V_fun_arg(3));
end

%% calculate second harmonic, AC E squared

% obtain averaged squared E field in the intrinsic
AC_Efield2 = AC_Efield .^ 2;
fit_AC_Efield2 = AC_Efield2(fit_t_index:end, :);

% average over absorber thickness, keep just time resolution

% to average the squared field is correct, even if it results in something
% that does not resemble a squared sin wave, for example one would expect
% this to touch zero but it does not if there is a phase distribution in
% the material
AC_Efield2_mean = trapz(absorber_x, AC_Efield2, 2) / absorber_thickness;
fit_AC_Efield2_mean = AC_Efield2_mean(fit_t_index:end);

% the phase dispersion is influencing the result, in order clarify its
% importance we can plot an hypotetical case where the electric field in
% every thin layer of perovskite is in phase, so that the sum is
% contructive

% this takes the amplitude of the squared AC electric field vs position
ac_amp_squared = trapz(s.t(fit_t_index:end), fit_AC_Efield .^ 2, 1) / (s.t(end) - s.t(fit_t_index));
% this averages the amplitudes obtained, on purpouse disregarding the
% importance of the phase distribution
coeff.ac_amp_squared_mean = trapz(absorber_x, ac_amp_squared) / absorber_thickness;

if demodulation
%    coeff.h2_mean = IS_EA_ana_demodulation(fit_t, fit_AC_Efield2_mean, s.par.V_fun_type, Vapp2_params);
    coeff.h2_mean = IS_EA_ana_demodulation(fit_t, fit_AC_Efield2_mean, s.par.V_fun_type, 2*s.par.V_fun_arg(3));
else
%    coeff.h2_mean = IS_EA_ana_fit(fit_t, fit_AC_Efield2_mean, s.par.E2_func);
    coeff.h2_mean = IS_EA_ana_fit(fit_t, fit_AC_Efield2_mean, s.par.V_fun_type, 2*s.par.V_fun_arg(3));
end

%% calculate ionic contribution

if s.par.mobseti && any(s.par.mucat) % if there was ion mobility, current due to ions have been calculated, fit it
    % calculate electric field due to ions in the oscillating solution
    tot_Efield_i_everywhere = dfana.calcFion(s);
    %s.par.e * cumtrapz(s.x(absorber_points), tot_i_matrix, 2) ./ s.par.dev.epp(absorber_points);
    tot_Efield_i = tot_Efield_i_everywhere(:, absorber_points);
    
    % the ionic electric field due to AC voltage does not oscillate around
    % the one in steady state solution (why?), so another way has to be used for
    % calculating baseline for ionic AC E field
    % the average field will be used as new DC field value
    DC_Efield_i = trapz(s.t(fit_t_index:end), tot_Efield_i(fit_t_index:end,:), 1) / (s.t(end) - s.t(fit_t_index));
    
    % calculate electric field due to AC variation of ions, first is
    % matrix, second is array
    AC_Efield_i = tot_Efield_i - DC_Efield_i;
    
    AC_ExDC_E_i = AC_Efield_i .* DC_Efield_i;
    fit_AC_ExDC_E_i = AC_ExDC_E_i(fit_t_index:end, :);

    AC_ExDC_E_i_mean = trapz(absorber_x, AC_ExDC_E_i, 2) / absorber_thickness;
    fit_AC_ExDC_E_i_mean = AC_ExDC_E_i_mean(fit_t_index:end);

    % remove bias and square the E field due to ions
    AC_Efield2_i = AC_Efield_i .^ 2;
    fit_AC_Efield2_i = AC_Efield2_i(fit_t_index:end, :);
    AC_Efield2_i_mean = trapz(absorber_x, AC_Efield2_i, 2) / absorber_thickness;
    fit_AC_Efield2_i_mean = AC_Efield2_i_mean(fit_t_index:end);

    % fit
    if demodulation
        coeff.i_h1_mean = IS_EA_ana_demodulation(fit_t, fit_AC_ExDC_E_i_mean, s.par.V_fun_type, s.par.V_fun_arg(3));
%        coeff.i_h2_mean = IS_EA_ana_demodulation(fit_t, fit_AC_Efield2_i_mean, s.par.V_fun_type, Vapp2_params);
        coeff.i_h2_mean = IS_EA_ana_demodulation(fit_t, fit_AC_Efield2_i_mean, s.par.V_fun_type, 2*s.par.V_fun_arg(3));
    else
        coeff.i_h1_mean = IS_EA_ana_fit(fit_t, fit_AC_ExDC_E_i_mean, s.par.V_fun_type, s.par.V_fun_arg(3));
 %       coeff.i_h2_mean = IS_EA_ana_fit(fit_t, fit_AC_Efield2_i_mean, s.par.E2_func);
        coeff.i_h2_mean = IS_EA_ana_fit(fit_t, fit_AC_Efield2_i_mean, s.par.V_fun_type, 2*s.par.V_fun_arg(3));
    end
end

if local_field
    
%% fit AC E field * DC E field in every space point, just if local_field is selected

    if demodulation
        % this would be much faster with vectorization of demodulation!!!
        for i = 1:absorber_points_num
            % if the profile is identically zero there's no need to run demodulation
            % this happens where DC component is zero, anyway this should
            % never happen
            if sum(abs(fit_AC_ExDC_E(:, i)))
                coeff.h1(:, i) = IS_EA_ana_demodulation(fit_t, fit_AC_ExDC_E(:, i), s.par.V_fun_type, s.par.V_fun_arg(3));
            else % if the profile is identically zero return amplitudes = zero and phase = not a number
                coeff.h1(:, i) = [0; 0; NaN];
            end
        end
    else
        % it's very unlikely to vectorize a fit :( and it's very slooow
        for i = 1:absorber_points_num
            % if the profile is identically zero the fit will return an error, anyway this should
            % never happen
            if sum(abs(fit_AC_ExDC_E(:, i)))
                coeff.h1(:, i) = IS_EA_ana_fit(fit_t, fit_AC_ExDC_E(:, i), s.par.V_fun_type, s.par.V_fun_arg(3));
            else % if the profile is identically zero return amplitudes = zero and phase = not a number
                coeff.h1(:, i) = [0; 0; NaN];
            end
        end
    end
    
%% fit AC E field in every space point and use it for first harmonic phase, just if local_field is selected     

    if demodulation
        % this would be much faster with vectorization of demodulation!!!
        for i = 1:absorber_points_num
            % if the profile is identically zero there's no need to run demodulation, anyway this should
            % never happen
            if sum(abs(fit_AC_Efield(:, i)))
                % notice the different Vapp_params, now for double frequency!
                coeff.ac(:, i) = IS_EA_ana_demodulation(fit_t, fit_AC_Efield(:, i), s.par.V_fun_type, s.par.V_fun_arg(3));
                % for obtaining the phase, both the AC profile or the
                % AC*DC profile can be used, but the former has less noise
                % this works unless DC component has negative sign
                coeff.h1(3, i) = coeff.ac(3, i);
                if DC_Efield(i) < 0 % not so common to happen
                    coeff.h1(2, i) = -coeff.h1(2, i);
                end
            else % if the profile is identically zero return amplitudes = zero and phase = not a number
                coeff.ac(:, i) = [0; 0; NaN];
            end
        end
    else
        % it's very unlikely to vectorize a fit :( and it's very slooow
        for i = 1:absorber_points_num
            % if the profile is identically zero the fit will return an error, anyway this should
            % never happen
            if sum(abs(fit_AC_Efield(:, i)))
                coeff.ac(:, i) = IS_EA_ana_fit(fit_t, fit_AC_Efield(:, i), s.par.V_fun_type, s.par.V_fun_arg(3));
                % for obtaining the phase, both the AC profile or the
                % AC*DC profile can be used, but the former has less noise
                % this works unless DC component has negative sign
                if DC_Efield(i) < 0
                    coeff.h1(3, i) = pi - coeff.ac(3, i);
                    coeff.h1(2, i) = -coeff.h1(2, i);
                else % shlould be quite always the case
                    coeff.h1(3, i) = coeff.ac(3, i);
                end
            else % if the profile is identically zero return amplitudes = zero and phase = not a number
                coeff.ac(:, i) = [0; 0; NaN];
            end
        end
    end
    
%% fit AC E field squared in every space point, just if local_field is selected

    if demodulation
        % this would be much faster with vectorization of demodulation!!!
        for i = 1:absorber_points_num
            % if the profile is identically zero there's no need to run demodulation, anyway this should
            % never happen
            if sum(abs(fit_AC_Efield2(:, i)))
                % notice the different Vapp_params, now for double frequency!
%                coeff.h2(:, i) = IS_EA_ana_demodulation(fit_t, fit_AC_Efield2(:, i), s.par.V_fun_type, Vapp2_params);
                coeff.h2(:, i) = IS_EA_ana_demodulation(fit_t, fit_AC_Efield2(:, i), s.par.V_fun_type, 2*s.par.V_fun_arg(3));
            else % if the profile is identically zero return amplitudes = zero and phase = not a number
                coeff.h2(:, i) = [0; 0; NaN];
            end
        end
    else
        % it's very unlikely to vectorize a fit :( and it's very slooow
        for i = 1:absorber_points_num
            % if the profile is identically zero the fit will return an error, anyway this should
            % never happen
            if sum(abs(fit_AC_Efield2(:, i)))
%                coeff.h2(:, i) = IS_EA_ana_fit(fit_t, fit_AC_Efield2(:, i), s.par.E2_func);
                coeff.h2(:, i) = IS_EA_ana_fit(fit_t, fit_AC_Efield2(:, i), s.par.V_fun_type, 2*s.par.V_fun_arg(3));
            else % if the profile is identically zero return amplitudes = zero and phase = not a number
                coeff.h2(:, i) = [0; 0; NaN];
            end
        end
    end

%% fit ionic field contributions in every space point, just if local_field is selected

    if demodulation
        % this would be much faster with vectorization of demodulation!!!
        for i = 1:absorber_points_num
            % if the profile is identically zero there's no need to run demodulation, anyway this should
            % never happen
            if sum(abs(fit_AC_ExDC_E_i(:, i)))
                coeff.i_h1(:, i) = IS_EA_ana_demodulation(fit_t, fit_AC_ExDC_E_i(:, i), s.par.V_fun_type, s.par.V_fun_arg(3));
            else % if the profile is identically zero return amplitudes = zero and phase = not a number
                coeff.i_h1(:, i) = [0; 0; NaN];
            end
            % if the profile is identically zero there's no need to run demodulation, anyway this should
            % never happen
            if sum(abs(fit_AC_Efield2_i(:, i)))
                % notice the different Vapp_params, now for double frequency!
%                coeff.i_h2(:, i) = IS_EA_ana_demodulation(fit_t, fit_AC_Efield2_i(:, i), s.par.V_fun_type, Vapp2_params);
                coeff.i_h2(:, i) = IS_EA_ana_demodulation(fit_t, fit_AC_Efield2_i(:, i), s.par.V_fun_type, 2*s.par.V_fun_arg(3));
            else % if the profile is identically zero return amplitudes = zero and phase = not a number
                coeff.i_h2(:, i) = [0; 0; NaN];
            end
        end
    else
        % it's very unlikely to vectorize a fit :( and it's very slooow
        for i = 1:absorber_points_num
            % if the profile is identically zero the fit will return an error, anyway this should
            % never happen
            if sum(abs(fit_AC_ExDC_E_i(:, i)))
                coeff.i_h1(:, i) = IS_EA_ana_fit(fit_t, fit_AC_ExDC_E_i(:, i), s.par.V_fun_type, s.par.V_fun_arg(3));
            else % if the profile is identically zero return amplitudes = zero and phase = not a number
                coeff.i_h1(:, i) = [0; 0; NaN];
            end
            % if the profile is identically zero the fit will return an error, anyway this should
            % never happen
            if sum(abs(fit_AC_Efield2_i(:, i)))
%                coeff.i_h2(:, i) = IS_EA_ana_fit(fit_t, fit_AC_Efield2_i(:, i), s.par.E2_func);
                coeff.i_h2(:, i) = IS_EA_ana_fit(fit_t, fit_AC_Efield2_i(:, i), s.par.V_fun_type, 2*s.par.V_fun_arg(3));
            else % if the profile is identically zero return amplitudes = zero and phase = not a number
                coeff.i_h2(:, i) = [0; 0; NaN];
            end
        end
    end
end

%% wrap phases to be nice to plot

coeff.h1_mean(3) = wrapToPi(coeff.h1_mean(3));
coeff.h1(3, :) = wrapToPi(coeff.h1(3, :));
coeff.ac(3, :) = wrapTo2Pi(coeff.ac(3, :));
coeff.i_h1_mean(3) = wrapTo2Pi(coeff.i_h1_mean(3));
coeff.i_h1(3, :) = wrapTo2Pi(coeff.i_h1(3, :));
coeff.h2_mean(3) = wrapToPi(coeff.h2_mean(3));
coeff.h2(3, :) = wrapToPi(coeff.h2(3, :));
coeff.i_h2_mean(3) = wrapToPi(coeff.i_h2_mean(3));
coeff.i_h2(3, :) = wrapToPi(coeff.i_h2(3, :));

%% plot solutions

if do_graphics % disable all this stuff if under parallelization or if explicitly asked to not plot any graphics
    Vapp_func = fun_gen(s.par.V_fun_type);
    Vapp = Vapp_func(s.par.V_fun_arg, s.t);
    x_nm = s.x(absorber_points) * 1e7;
   
    % fourth value of Vapp_params is pulsatance
    figure('Name', ['Single EA, Int ' num2str(s.par.int1)...
        ' Freq ' num2str(s.par.V_fun_arg(3)) ' - 1st harmonic averaged over position'], 'NumberTitle', 'off');
        yyaxis right
        hold off
        h(1) = plot(s.t, Vapp, 'r');
        legend_array = "Applied Voltage";
        ylabel('Applied voltage [V]');

        yyaxis left
        hold off
        h(2) = plot(s.t, AC_ExDC_E_mean, 'b');
        hold on
        h(3) = plot(fit_t, Vapp_func([coeff.h1_mean(1),coeff.h1_mean(2),s.par.V_fun_arg(3),coeff.h1_mean(3)], fit_t), 'kx-');
        legend_array = [legend_array, "Averaged E field", "E field fit"];
%         if s.par.mucat % if there was ion mobility, E field due to ions have been calculated, plot stuff
%             h(4) = plot(s.t, AC_ExDC_E_i_mean, 'g-');
%             h(5) = plot(fit_t, Vapp_func([coeff.i_coeff_1h_mean(1),coeff.i_coeff_1h_mean(2),s.par.V_fun_arg(3),coeff.i_coeff_1h_mean(3)], fit_t), 'kx-');
%             legend_array = [legend_array, "Ionic E field", "Ionic E field fit"];
%         end
        ylabel('E_{AC} x E_{DC} amplitude [V^2/cm^2]');
        xlabel('Time [s]');
        legend(h, legend_array);   
        legend boxoff
        clear h
        
    figure('Name', ['Single EA, Int ' num2str(s.par.int1)...
        ' Freq ' num2str(s.par.V_fun_arg(3)) ' - 2nd harmonic averaged over position'], 'NumberTitle', 'off');
        yyaxis right
        hold off
        h(1) = plot(s.t, Vapp, 'r');
        legend_array = "Applied Voltage";
        ylabel('Applied voltage [V]');

        yyaxis left
        hold off
        h(2) = plot(s.t, AC_Efield2_mean, 'b');
        hold on
        h(3) = plot(fit_t, Vapp_func([coeff.h2_mean(1),coeff.h2_mean(2),2*s.par.V_fun_arg(3),coeff.h2_mean(3)], fit_t), 'kx-');

        legend_array = [legend_array, "Averaged E^2 field", "E^2 field fit"];
%         if s.par.mucat % if there was ion mobility, E2 field due to ions have been calculated, plot stuff
%             h(4) = plot(s.t, AC_Efield2_i_mean, 'g-');
%             h(5) = plot(fit_t, s.par.E2_func(i_coeff_2h_mean, fit_t), 'kx-');
%             legend_array = [legend_array, "Ionic E^2 field", "Ionic E^2 field fit"];
%         end
        ylabel('E_{AC}^2 amplitude [V^2/cm^2]');
        xlabel('Time [s]');
        legend(h, legend_array);
        legend boxoff
        clear h
        
        

    if local_field

        figure('Name', ['Single EA, Int ' num2str(s.par.int1)...
            ' Freq ' num2str(s.par.V_fun_arg(3)) ' - 1st harmonic vs position'], 'NumberTitle', 'off');
            yyaxis left
            hold off
            i = 1;
            h(i) = plot(x_nm, rad2deg(coeff.h1(3, :)), 'b');
            hold on
            legend_array = "Phase";
%             if s.par.mucat % if there was ion mobility
%                 i = i + 1; h(i) = plot(x_nm, rad2deg(i_coeff_1h(3, :)), 'b--');
%                 legend_array = [legend_array, "Ionic phase"];
%             end
            ylabel('E_{AC} x E_{DC} phase [deg]');
            
            yyaxis right
            hold off
            i = i + 1; h(i) = plot(x_nm, coeff.h1(2, :), 'r');
            hold on
            legend_array = [legend_array, "Amplitude"];
%             if s.par.mucat % if there was ion mobility
%                 i = i + 1; h(i) = plot(x_nm, i_coeff_1h(2, :), 'r--');
%                 legend_array = [legend_array, "Ionic amplitude"];
%             end
            ylabel('E_{AC} x E_{DC} amplitude [V^2/cm^2]');
            xlabel('Position [nm]');
            legend(h, legend_array);
            legend boxoff
            clear h

        fig_2H_vs_x = figure('Name', ['Single EA, Int ' num2str(s.par.int1)...
            ' Freq ' num2str(s.par.V_fun_arg(3)) ' - 2nd harmonic vs position'], 'NumberTitle', 'off');
            yyaxis left
            hold off
            i = 1;
            h(i) = plot(x_nm, rad2deg(coeff.h2(3, :)), 'b');
            hold on
            legend_array = "Phase";
%             if s.par.mucat % if there was ion mobility
%                 i = i + 1; h(i) = plot(x_nm, rad2deg(i_coeff_2h(3, :)), 'b--');
%                 legend_array = [legend_array, "Ionic phase"];
%             end
            ylabel('E_{AC}^2 phase [deg]');
            
            yyaxis right
            hold off
            i = i + 1; h(i) = plot(x_nm, coeff.h2(2, :), 'r');
            hold on
            legend_array = [legend_array, "Amplitude"];
%             if s.par.mucat % if there was ion mobility
%                 i = i + 1; h(i) = plot(x_nm, i_coeff_2h(2, :), 'r--');
%                 legend_array = [legend_array, "Ionic amplitude"];
%             end
            ylabel('E_{AC}^2 amplitude [V^2/cm^2]');
            xlabel('Position [nm]');
            legend(h, legend_array);
            legend boxoff
            clear h

        figure('Name', ['Single EA, Int ' num2str(s.par.int1)...
            ' Freq ' num2str(s.par.V_fun_arg(3)) ' - AC and DC fields intensity vs position'], 'NumberTitle', 'off');
            yyaxis left
            hold off
            i = 1;
            h(i) = plot(x_nm, DC_Efield, 'b');
            hold on
            legend_array = "DC E field intensity";
            ylabel('E_{DC} intensity [V/cm]');
            
            yyaxis right
            hold off
            i = i + 1; h(i) = plot(x_nm, coeff.ac(2, :), 'r');
            hold on
            legend_array = [legend_array, "AC E field amplitude"];
            ylabel('E_{AC} amplitude [V/cm]');
            xlabel('Position [nm]');
            legend(h, legend_array);
            legend boxoff
            clear h
    end
    if ~ismissing(savefig_dir)
        if 7~=exist(savefig_dir,'dir')
            mkdir(savefig_dir)
        end
        pathname = [char(savefig_dir) filesep char(inputname(1)) '-2H'];
        savefig(fig_2H_vs_x, pathname, 'compact')
        saveas(fig_2H_vs_x, [pathname '.png'])
    end
end

%------------- END OF CODE --------------
