function [extracting_array, index_early] = CE_ISstep_single_analysis(asymstruct_CE_ISstep, Int, subtract_baseline, do_graphics)
%CE_ISSTEP_SINGLE_ANALYSIS - Calculate charge extracted by Charge Extraction (CE)
% and Impedance Spectroscopy approximated by a step (ISstep) integrating the current
%
% Syntax:  [extracting_array, index_early] = CE_ISstep_single_analysis(asymstruct_CE_ISstep, Int, subtract_baseline)
%
% Inputs:
%   ASYMSTRUCT_CE_ISSTEP - a struct with a solution being perturbed by CE
%     or ISstep
%   INT - original light intensity is needed in input as the solution provided
%     in input for CE is in dark (this is not true for ISstep)
%   SUBTRACT_BASELINE - logical indicating if the residual current should
%     be subtracted, usually true for ISstep and false for CE
%   DO_GRAPHICS - logical, whether to calculate the ionic contribution and graph the solution
%
% Outputs:
%   EXTRACTING_ARRAY - array of cumulative charge extracted at each time point;
%   INDEX_EARLY - index of the previous array indicating the end of the early extraction, time at which the electronic charge have been extracted but the ionic defects did not move yet.
%
% Example:
%   CE_ISstep_single_analysis(CE_single_exec(asymmetricize(ssol_i_light, 1), 1e-3, 1), 1, false, true)
%     integrate the current and plot a CE solution
%   CE_ISstep_single_analysis(ISstep_single_exec(asymmetricize(ssol_i_light, 1), 1e-2, 1, 1e-3, 200), 1, true, true)
%     integrate the current and plot a ISstep solution
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also CE_single_exec, CE_full_exec, ISstep_single_exec, ISstep_full_exec.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

% increase graphics font size
set(0, 'defaultAxesFontSize', 24);
% set image dimension
set(0, 'defaultfigureposition', [0, 0, 1000, 750]);
% set line thickness
set(0, 'defaultLineLineWidth', 2);

% evil shortcut
s = asymstruct_CE_ISstep;

% get the current in Ampere
Jamp = -s.Jn ./ 1000;

%% do the integral over time
if subtract_baseline
    % if the Vapp at the end is kept constant, like in CE and ISstep,
    % likely the final value (average of last 3 points) of the current has to be subtracted
    extracting_array = cumtrapz(s.t, Jamp - mean(Jamp(end-2:end)));
else
    % if the Vapp is changing in some other way, maybe the current does not
    % stabilize, so can't be subtracted
    extracting_array = cumtrapz(s.t, Jamp);
end
extracting_charges = extracting_array(end); % final value of the integral

%% get just the early extracted charges
% I hope that this will never match a point before the current peak
[deltaJ, peak_index] = min(s.Jn);
% subsetting just current after peak for having the possibility to get the
% time of the decay
J_after_peak = transpose(s.Jn(peak_index:end));
t_after_peak = s.t(peak_index:end);
[~, index_tenth_after] = min(abs(J_after_peak - deltaJ / 10)); % index of very low current after decay
t_tenth = t_after_peak(index_tenth_after);
t_early_approx = t_tenth * 5; % early time extraction ends some times after the end of current decay
[~, index_early] = min(abs(s.t - t_early_approx)); % get the index of an actual point close to that time
t_early = s.t(index_early); % get an actual time point where the charge is stabilized after fast decay

disp([mfilename ' - Extracting charges ' num2str(extracting_charges)...
    ' C/cm2; Early time charges ' num2str(extracting_array(index_early)) ' C/cm2'])

if do_graphics

%% calculate ionic contribution
    % calculate displacement current due to ionic movement
    % based on the same function from ISwave_single_analysis
    if s.p.mui % if there was ion mobility, current due to ions have been calculated, fit it
        % Intrinsic points logical array
        itype_points= (s.x >= s.p.tp & s.x <= s.p.tp + s.p.ti);
        % subtract background ion concentration, for having less noise in trapz
        i_matrix = s.sol(:, :, 3) - s.p.NI;
        % calculate electric field due to ions
        Efield_i = s.p.e * cumtrapz(s.x, i_matrix, 2) / s.p.eppi;
        % an average would be enough if the spatial mesh was homogeneous in the
        % intrinsic, indeed I have to use trapz for considering the spatial mesh
        Efield_i_mean = trapz(s.x(itype_points), Efield_i(:, itype_points), 2) / s.p.ti;
        % calculate displacement current due to ions
        Ji_disp = s.p.eppi * gradient(Efield_i_mean, s.t); % in Amperes
        Ji_disp = 1000 * Ji_disp; % from A to mA
    else
        Ji_disp = NaN;
    end
    
    [~, ~, U] = pinana(asymstruct_CE_ISstep);
    
%% plot solutions
    figure('Name', ['Single CE or IS at light intensity ' num2str(Int)], 'NumberTitle', 'off');
        yyaxis right
        hold off
        h(1) = plot(s.t, -extracting_array);
        legend_text = "Integrated Charge";
        ylabel('Charge [C/cm^2]');
        ax = gca;
        ax.XScale = 'log'; % for putting the scale in log

        yyaxis left
        hold off
        h(2) = plot(s.t, s.Jn, 'b--');
        hold on
        legend_text = [legend_text, "Total Current"];
        if asymstruct_CE_ISstep.p.mui
            h(3) = plot(s.t, -Ji_disp, 'b-');
            current_ylim = get(gca, 'ylim');
            plot([t_early, t_early], current_ylim, 'k--'); % draw vertical line for time of extraction of early current
            legend_text = [legend_text, "Ionic Displ. Curr."];
        end
        ylabel('Current [mA/cm^2]');
        hold off
        xlabel('Time [s]');
        legend(h, legend_text);
        hold off
        
    figure('Name', ['Single CE or IS at light intensity ' num2str(Int)], 'NumberTitle', 'off');
        if asymstruct_CE_ISstep.p.mui
            yyaxis right
            hold off
            h(1) = plot(s.t, -Ji_disp, 'b-');
            legend_text = "Ionic Displ. Curr.";
            ylabel('Ionic Current [mA/cm^2]');
            ax = gca;
            ax.XScale = 'log'; % for putting the scale in log
        end

        yyaxis left
        hold off
        h(2) = plot(s.t, s.Jn, 'b--');
        hold on
        legend_text = [legend_text, "Total Current"];
        h(3) = plot(s.t, U-min(U), 'k:');
        legend_text = [legend_text, "Recomb. Curr. - baseline"];

        ylabel('Current [mA/cm^2]');
        hold off
        xlabel('Time [s]');
        legend(h, legend_text);
        hold off
end

%------------- END OF CODE --------------
