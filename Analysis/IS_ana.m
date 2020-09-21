function coeff = IS_ana(s, minimal_mode, demodulation)
%IS_ANA - Calculate impedance (reactance and resistance) and phase by Impedance Spectroscopy (IS) with oscillating voltage
%
% Syntax:  IS_ana(s, minimal_mode, demodulation)
%
% Inputs:
%   S - a struct with a solution being perturbed by an
%     oscillating voltage, as generated from doIS_EA
%   MINIMAL_MODE - logical, when true graphics does not get created and
%     IS_subtracting_analysis does not get launched, useful when
%     launched under parallelization
%   DEMODULATION - logical, get phase via demodulation instead of using a fitting
%
% Outputs:
%   coeff.Jtot - array of: background total current, half peak to peak amplitude of
%     oscillation and phase of total electronic current
%   coeff.ion_disp - array of: background ionic current, half peak to peak amplitude of
%     oscillation and phase of ionic displacement current
%   coeff.r - array of: background recombination current, half peak to peak amplitude of
%     oscillation and phase of recombinating charge per unit time
%   coeff.np_dt - array of: background accumulating current, half peak to peak amplitude of
%     oscillation and phase of accumulating current, this is the real
%     capacitive current, obtained comparing the free charges profiles at
%     different times
%  
% Example:
%   IS_ana(ssol_i_1S_SR_is_100mHz_2mV, false, true)
%     plot current profile, reference profiles and calculate the phase using demodulation approach
%
% Other m-files required: IS_ana_subtracting, IS_EA_ana_fit, IS_EA_ana_demodulation, dfana
% Subfunctions: none
% MAT-files required: none
%
% See also IS_script, doIS_EA, IS_ana_subtracting, IS_EA_ana_demodulation, IS_EA_ana_fit.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: May 2020

%------------- BEGIN CODE --------------

% verify if the simulation broke, in that case return just NaNs
if size(s.u, 1) < s.par.tpoints
    coeff.Jtot = [NaN, NaN, NaN];
    coeff.ion_disp = [NaN, NaN, NaN];
    coeff.cat_disp = [NaN, NaN, NaN];
    coeff.ani_disp = [NaN, NaN, NaN];
    coeff.r = [NaN, NaN, NaN];
    coeff.np_dt = [NaN, NaN, NaN];
%     c_n_noionic = [NaN, NaN, NaN];
    return;
end

%% get current profiles

% s.par.V_fun_arg(3) is frequency
% round _should_ not be needed here
periods = round(s.par.tmax * s.par.V_fun_arg(3));

% here is critical that exactely an entire number of periods is provided to
% IS_single_demodulation
fit_t_index = round((s.par.tpoints - 1) * floor(periods / 2) / periods) + 1;
fit_t = s.t(fit_t_index:end)';

% current profile to be analyzed
J = dfana.calcJ(s);
Jtotx = J.tot(:, end);
fit_J = Jtotx(fit_t_index:end); % in Ampere

% recombination flux
r = dfana.calcr(s);
r_time_profile = trapz(s.x, r.tot, 2) * s.par.e;
fit_r = r_time_profile(fit_t_index:end); % in Ampere

% accumulating current profile to be analyzed
[np_dt, ~] = IS_ana_subtracting(s);
fit_np_dt = np_dt(fit_t_index:end); % in Ampere

%% remove some tilting from fit_J

% % to get better fit and better demodulation in case of
% % unstabilized solutions (not usually the case). In case of noisy solutions this could work badly
% % but should not affect too much the fitting/demodulation
% delta_t = fit_t(end) - fit_t(1);
% % this assumes that the first and last point are at the same point in the oscillating voltage
% delta_J = fit_J(end) - fit_J(1);
% tilting = delta_J/delta_t;
% t_middle = fit_t(round(end/2));
% % because of this, the bias value that will be obtained from the fit/demodulation is not going to be correct
% fit_J_flat = fit_J - tilting * (fit_t - t_middle);

%% extract parameters from current profiles

% s.par.V_fun_arg(3) is frequency
if demodulation
    coeff.Jtot = IS_EA_ana_demodulation(fit_t, fit_J, s.par.V_fun_type, s.par.V_fun_arg(3));
    coeff.r = IS_EA_ana_demodulation(fit_t, fit_r, s.par.V_fun_type, s.par.V_fun_arg(3));
    coeff.np_dt = IS_EA_ana_demodulation(fit_t, fit_np_dt', s.par.V_fun_type, s.par.V_fun_arg(3));
else
    coeff.Jtot = IS_EA_ana_fit(fit_t, fit_J, s.par.V_fun_type, s.par.V_fun_arg(3));
    coeff.r = IS_EA_ana_fit(fit_t, fit_r, s.par.V_fun_type, s.par.V_fun_arg(3));
    coeff.np_dt = IS_EA_ana_fit(fit_t, fit_np_dt, s.par.V_fun_type, s.par.V_fun_arg(3));
end

%% calculate ionic contribution

if s.par.mobseti && any(s.par.mucat) % if there was ion mobility, current due to ions have been calculated, fit it
    % Intrinsic points logical array
    absorber_points = logical(s.par.dev.g0);
%     % subtract background ion concentration, for having less noise in trapz
%     i_matrix = s.u(:, :, 4) - s.par.dev.Ncat;
%     % get absorber thickness
     absorber_x = s.x(absorber_points);
     absorber_thickness = (absorber_x(end)-absorber_x(1));
%     % calculate electric field due to ions
%     Efield_i = -s.par.e * cumtrapz(absorber_x, i_matrix(:,absorber_points), 2) ./ s.par.dev.epp(absorber_points);
%     % an average would be enough if the spatial mesh was homogeneous in the
%     % intrinsic, indeed I have to use trapz for considering the spatial mesh
%     Efield_i_per_epp = Efield_i .* s.par.dev.epp(absorber_points);
%     Efield_i_per_epp_mean = trapz(s.x(absorber_points), Efield_i_per_epp, 2) / absorber_thickness;
%     % calculate displacement current due to ions
%     Ji_disp = gradient(Efield_i_per_epp_mean, s.t); % in Amperes
%    fit_Ji = Ji_disp(fit_t_index:end); % in Ampere

    F_ion = dfana.calcFion(s);
    F_ion_per_epp = -s.par.e .* F_ion(:,absorber_points) .* s.par.dev.epp(absorber_points) .* s.par.epp0;
    F_ion_per_epp_mean = trapz(s.x(absorber_points), F_ion_per_epp, 2) / absorber_thickness;
    ion_disp = gradient(F_ion_per_epp_mean, s.t); % in Amperes
    fit_ion_disp = ion_disp(fit_t_index:end); % in Ampere

    if demodulation
%         coeff.ion_disp = IS_EA_ana_demodulation(fit_t, fit_Ji, s.par.V_fun_type, s.par.V_fun_arg(3));
        coeff.ion_disp = IS_EA_ana_demodulation(fit_t, fit_ion_disp, s.par.V_fun_type, s.par.V_fun_arg(3));
    else
%         coeff.ion_disp = IS_EA_ana_fit(fit_t, fit_Ji, s.par.V_fun_type, s.par.V_fun_arg(3));
        coeff.ion_disp = IS_EA_ana_fit(fit_t, fit_ion_disp, s.par.V_fun_type, s.par.V_fun_arg(3));
    end
    
        
    if size(s.u,3) > 4
        [~,~,x,par,dev,~,~,a,c,~] = dfana.splitsol(s);

        rho_cat = c - dev.Ncat;
        rho_ani = - a + dev.Nani;
        
        F_cat = cumtrapz(x, rho_cat, 2)./(dev.epp*par.epp0);
        F_ani = cumtrapz(x, rho_ani, 2)./(dev.epp*par.epp0);
        
        F_cat_per_epp = -s.par.e .* F_cat(:,absorber_points) .* s.par.dev.epp(absorber_points) .* s.par.epp0;
        F_cat_per_epp_mean = trapz(s.x(absorber_points), F_cat_per_epp, 2) / absorber_thickness;
        cat_disp = gradient(F_cat_per_epp_mean, s.t); % in Amperes
        fit_cat_disp = cat_disp(fit_t_index:end); % in Ampere
    
        F_ani_per_epp = -s.par.e .* F_ani(:,absorber_points) .* s.par.dev.epp(absorber_points) .* s.par.epp0;
        F_ani_per_epp_mean = trapz(s.x(absorber_points), F_ani_per_epp, 2) / absorber_thickness;
        ani_disp = gradient(F_ani_per_epp_mean, s.t); % in Amperes
        fit_ani_disp = ani_disp(fit_t_index:end); % in Ampere
        if demodulation
            %         coeff.ion_disp = IS_EA_ana_demodulation(fit_t, fit_Ji, s.par.V_fun_type, s.par.V_fun_arg(3));
            coeff.cat_disp = IS_EA_ana_demodulation(fit_t, fit_cat_disp, s.par.V_fun_type, s.par.V_fun_arg(3));
            coeff.ani_disp = IS_EA_ana_demodulation(fit_t, fit_ani_disp, s.par.V_fun_type, s.par.V_fun_arg(3));
        else
            %         coeff.ion_disp = IS_EA_ana_fit(fit_t, fit_Ji, s.par.V_fun_type, s.par.V_fun_arg(3));
            coeff.cat_disp = IS_EA_ana_fit(fit_t, fit_cat_disp, s.par.V_fun_type, s.par.V_fun_arg(3));
            coeff.ani_disp = IS_EA_ana_fit(fit_t, fit_ani_disp, s.par.V_fun_type, s.par.V_fun_arg(3));
        end
    else
        coeff.cat_disp = coeff.ion_disp;
    end
    
%     
%     % calculate electronic current subtracting ionic contribution
%     Jn_noionic = Jtotx - Ji_disp; % in Ampere

else % if no ionic mobility is present, report NaNs
    coeff.ion_disp = [NaN, NaN, NaN];
%     Jn_noionic = NaN;
end

%% plot solutions

if ~minimal_mode % disable all this stuff if under parallelization or if explicitly asked to not plot any graphics

    Vapp_func = fun_gen(s.par.V_fun_type);
    func_fixedfreq = @(coeff,t) Vapp_func([coeff(1), coeff(2), s.par.V_fun_arg(3), coeff(3)], t);
    
    % in phase electronic current
    Jn_inphase = func_fixedfreq([coeff.Jtot(1), coeff.Jtot(2)*cos(coeff.Jtot(3)), 0], s.t);
    % out of phase electronic current
    Jn_quadrature = func_fixedfreq([coeff.Jtot(1), coeff.Jtot(2)*sin(coeff.Jtot(3)), pi/2], s.t);

%     if s.par.mobseti && any(s.par.mucat)
%         fit_Jn_noionic = Jn_noionic(fit_t_index:end);
%         if demodulation
%             c_n_noionic = IS_EA_ana_demodulation(fit_t, fit_Jn_noionic, s.par.V_fun_type, s.par.V_fun_arg(3));
%         else
%             c_n_noionic = IS_EA_ana_fit(fit_t, fit_Jn_noionic, s.par.V_fun_type, s.par.V_fun_arg(3));
%         end
% 
% %         % in phase electronic current
% %         Jn_noionic_inphase = func_fixedfreq([n_noionic_coeff(1), n_noionic_coeff(2)*cos(n_noionic_coeff(3)), 0], s.t);
% %         % out of phase electronic current
% %         Jn_noionic_quadrature = func_fixedfreq([n_noionic_coeff(1), n_noionic_coeff(2)*sin(n_noionic_coeff(3)), pi/2], s.t);
% 
%     else
% %         Jn_noionic_inphase = NaN;
% %         Jn_noionic_quadrature = NaN;
%         c_n_noionic = [NaN, NaN, NaN];
%     end

    % fourth value of Vapp_params is pulsatance
    figure('Name', ['Single IS, Int ' num2str(s.par.int1) ' Freq '...
        num2str(s.par.V_fun_arg(3))], 'NumberTitle', 'off');
        yyaxis right
        hold off
        i=1;
        h(i) = plot(s.t, Vapp_func(s.par.V_fun_arg, s.t), 'r', 'LineWidth', 2);
        legend_array = "Applied Voltage";
        ylabel('Applied voltage [V]');

        yyaxis left
        hold off
        i=i+1; h(i) = plot(s.t, Jtotx * 1000, 'k-', 'LineWidth', 2); % mA
        legend_array = [legend_array, "Current"];
        hold on
        i=i+1; h(i) = plot(s.t, r_time_profile * 1000, 'k--'); % mA
        legend_array = [legend_array, "Recombination current"];
        i=i+1; h(i) = plot(s.t, np_dt * 1000, 'b:', 'LineWidth', 2); % mA
        legend_array = [legend_array, "Accumulating current"];
        i=i+1; h(i) = plot(fit_t, func_fixedfreq(coeff.Jtot, fit_t) * 1000, 'kx-'); % mA
        legend_array = [legend_array, "Fit of Current"];
        i=i+1; h(i) = plot(s.t, Jn_inphase*1000, 'm-', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 7); % mA
        legend_array = [legend_array, "In phase J"];
        i=i+1; h(i) = plot(s.t, Jn_quadrature*1000, 'm-', 'LineWidth', 1, 'Marker', 'x', 'MarkerSize', 7); % mA
        legend_array = [legend_array, "Out of phase J"];
        if s.par.mobseti && any(s.par.mucat) % if there was ion mobility, current due to ions have been calculated, plot stuff
            i=i+1; h(i) = plot(s.t, ion_disp * 1000, 'g--', 'LineWidth', 2); % mA
            legend_array = [legend_array, "Ionic displacement current"];
        end
        ylabel('Current [mA/cm^2]');
        hold off
        xlabel('Time [s]');
        legend(h, legend_array);
        hold off
end

%------------- END OF CODE --------------
