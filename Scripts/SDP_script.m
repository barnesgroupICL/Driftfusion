function sdpsol = SDP_script(sol_ini, tdwell_arr, Vjump, bias_source, bias_int, dwell_mode, pulse_source, pulse_int, pulse_tmax, pulse_mobile_ions)
%SDP_SCRIPT - an alternative script for Step-Dwell-Probe simulations.
%
% Syntax:  sdpsol = SDP_script(sol_ini, tdwell_arr, Vjump, bias_source, bias_int, pulse_source, pulse_int, pulse_tmax, pulse_mobile_ions)
%
% Inputs:
%   SOL_INI - a single struct created by DF.
%   TDWELL_ARR - an array of floats indicating the dwell time for each SDP
%     simulation
%   VJUMP - a float, the amplitude of the voltage jump to be performed
%     before the dwell step
%   BIAS_SOURCE - an integer, which illumination source should be turn on 
%     before the dwell step, as defined in Optical/beerlambert.m
%   BIAS_INT - a float, the intensity of the bias illumination. For the AM15
%     source the 1 is equivalent to 1 sun illumination
%   DWELL_MODE - a string, when 'separated' each dwell step is calculated
%     from time zero up to the requested time; when 'sequential' a faster
%     method is used in which each dwell is an extension of the previous
%     simulation; when 'unique' all the dwell solutions are extracted from
%     an unique simulation which includes all the requested dwell times.
%     Sometimes the three methods could give slightly different results. In
%     these cases increase the number of points in the time mesh, editing
%     the code where par_*.tpoints is set.
%   PULSE_SOURCE - an integer, which illumination source should be turn on 
%     after the dwell step (the pulse, a.k.a. probe), as defined in Optical/beerlambert.m
%   PULSE_INT - a float, the intensity of the pulse (a.k.a. probe) illumination. For the AM15
%     source the 1 is equivalent to 1 sun illumination
%   PULSE_TMAX - a float, the duration of the pulse (probe) step at the end
%     of which the SDP numeric output is evaluated. Experimentally this
%     number should be around one microsecond, so that the electronic
%     charge had time to adapt at the pulse conditions but the ionic
%     charges are still at the status caused by the dwell step.
%     A large value can be set if the PULSE_MOBILE_IONS option is set to
%     false.
%   PULSE_MOBILE_IONS - a boolean, whether the ions should be mobile or not
%     in the pulse (probe) step. As in this step we're usually interested in
%     the electronic behaviour at times when the ions still didn't move, we
%     can as well freeze the ions.
%
% Outputs:
%   SDPSOL - a structure, including: Jdk, the current value at the end of the dwell
%     step; Jtr, the current profile during the pulse step minus the
%     Jdk current; t_Jtr, the time axis during the pulse step; and a few
%     parameters as supplied in the input arguments.
%
% Example:
%   sdpsol_spiro_01sun = SDP_script(soleq_spiro.ion, [logspace(-8,3,56)], 0.6, 1, 0.1, 2, 5.12, 1e-3, true);
%     perform 56 SDP simulations using dwell times from 10 ns to 1000 s,
%     using a voltage step of 0.6 V, illuminating from the primary light
%     source (AM15 by default) with intensity of 0.1 sun during the dwell
%     step, adding a secondary illumination (a red laser by default) with
%     an intensity of 5.12 as the pulse with a duration of 1 ms, which is
%     too long for ensuring that the ions didn't move, so that the first
%     "true" value instructs the solver to freeze the ions during this
%     step.
%
% Other m-files required: df, dfana, runDfFourTrunks
% Subfunctions: runDfLoosenTolAndReduceMaxStepFactor, calcJtrJdk,
%   lightRampUp, lightConstant
% MAT-files required: none
%
% See also df, anasdp, SDP_list_plot, SDP_script_exporter, doSDP, runDfFourTrunks.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------
disp('Starting SDP protocol')

jump_time = 1e-10;

% Get parameters
par_jump = sol_ini.par;

% tmax practically sets the minimum dwell time
par_jump.tmesh_type = 1;
par_jump.tmax = jump_time;
par_jump.tpoints = 20;
% Voltage jump
par_jump.V_fun_type = 'sweep';
% sol_ini.par.Vapp usually is zero
par_jump.V_fun_arg = [sol_ini.par.Vapp, Vjump, par_jump.tmax];
% freeze the ions
par_jump.mobseti = 0;

if bias_int
    % the initial intensity is zero as we assume that the starting solution
    % is in dark
    par_jump = lightRampUp(par_jump, bias_source, 0, bias_int);
end

disp('Initial voltage (and illumination) jump')
sol_jump = runDfLoosenTolAndReduceMaxStepFactor(sol_ini, par_jump);

par_dwell = lightConstant(par_jump, bias_source, bias_int);
% unset Vjump parameters
par_dwell.V_fun_type = 'constant';
par_dwell.V_fun_arg = Vjump;
% thaw ions
par_dwell.mobseti = 1;
par_dwell.tmesh_type = 2;
% unexpectedly, the number of tpoints in the dwell stage is important
par_dwell.tpoints = 50;

% usually, two different light sources are used for the bias and the probe,
% but if this is not the case, the ramp up should start from the already
% applied value
if pulse_source == bias_source
    pre_pulse_int = bias_int;
else
    pre_pulse_int = 0;
end
par_pulseTurnOn = lightRampUp(par_dwell, pulse_source, pre_pulse_int, pulse_int);
par_pulseTurnOn.mobseti = 0;
par_pulseTurnOn.tmesh_type = 1;
par_pulseTurnOn.tmax = jump_time;

par_pulseKeepOn = lightConstant(par_pulseTurnOn, pulse_source, pulse_int);
% we should not use pulses long enough to have ions moving, so we could as
% well freeze them during the probe pulse, for helping the solver
par_pulseKeepOn.mobseti = pulse_mobile_ions;
par_pulseKeepOn.tmesh_type = 2;
par_pulseKeepOn.t0 = jump_time;
par_pulseKeepOn.tmax = pulse_tmax;
par_pulseKeepOn.tpoints = 100;

% Preallocate memory
Jtr = zeros(par_pulseKeepOn.tpoints, length(tdwell_arr));
t_Jtr = zeros(par_pulseKeepOn.tpoints, length(tdwell_arr));
Jdk = zeros(1,length(tdwell_arr));

switch dwell_mode
    case 'separated'
        unique_dwell = false;
        sequential_dwell = false;
    case 'sequential'
        unique_dwell = false;
        sequential_dwell = true;
    case 'unique'
        unique_dwell = true;
end

if unique_dwell
    par_dwell.tmax = max(tdwell_arr);
    par_dwell.t0 = min(tdwell_arr)/10;
    % unexpectedly, the number of tpoints in the dwell stage is important
    par_dwell.tpoints = 50;
    % generate tmesh
    tmesh = meshgen_t(par_dwell);
    % custom time mesh
    par_dwell.tmesh_type = 5;
    % add desired points to tmesh
    par_dwell.t = unique([tdwell_arr, tmesh]);
    par_dwell.tpoints = length(par_dwell.t);
    sol_dwell_unique = df(sol_jump, par_dwell);
    sol_dwell.x = sol_dwell_unique.x;
else
    firstRun = true;
    brokenRun = true;
end

for i = 1:length(tdwell_arr)
    if unique_dwell
        tpoints_before_requested_time = sol_dwell_unique.t <= tdwell_arr(i);
        sol_dwell.t = sol_dwell_unique.t(tpoints_before_requested_time);
        sol_dwell.u = sol_dwell_unique.u(tpoints_before_requested_time,:,:);
        sol_dwell.par = sol_dwell_unique.par;
        sol_dwell.par.tpoints = length(sol_dwell.t);
    else
        disp(['SDP tdwell = ' num2str(tdwell_arr(i)) ' - Dwell'])
        if ~firstRun && sequential_dwell
            par_dwell.tmax = tdwell_arr(i) - tdwell_arr(i-1);
            par_dwell.t0 = par_dwell.tmax / 20;
            sol_dwell = df(sol_dwell, par_dwell);
            if size(sol_dwell.u,1) == par_dwell.tpoints
                brokenRun = false;
            end
        end
        if firstRun || brokenRun || ~sequential_dwell
            par_dwell.tmax = tdwell_arr(i);
            % we care just about the final timepoint of the dwell, so a too
            % large t0 should not be a problem here
            par_dwell.t0 = par_dwell.tmax / 20;
            sol_dwell = df(sol_jump, par_dwell);
            firstRun = false;
            if size(sol_dwell.u,1) == par_dwell.tpoints
                brokenRun = false;
            end
        end
        if brokenRun
            sol_dwell = runDfFourTrunks(sol_jump, par_dwell);
        end
    end

    
    disp(['SDP tdwell = ' num2str(tdwell_arr(i)) ' - Pulse turn on'])
    sol_pulseTurnOn = runDfLoosenTolAndReduceMaxStepFactor(sol_dwell, par_pulseTurnOn);
    
    disp(['SDP tdwell = ' num2str(tdwell_arr(i)) ' - Pulse keep on'])
    sol_pulseKeepOn = df(sol_pulseTurnOn, par_pulseKeepOn);
    if size(sol_pulseKeepOn.u,1) ~= par_pulseKeepOn.tpoints
        sol_pulseKeepOn = runDfFourTrunks(sol_pulseTurnOn, par_pulseKeepOn);
    end
    [Jtr_temp, Jdk(i)] = calcJtrJdk(sol_dwell, sol_pulseKeepOn);
    t_Jtr_temp = sol_pulseKeepOn.t';
    % runDfFourTrunks can output a solution with more time points than
    % requested, but anasdp works only on matrices, not on cells with
    % arrays of different lengths, so the interpolated output is used
    if all(size(Jtr_temp) == size(Jtr(:,i)))
        Jtr(:,i) = Jtr_temp;
        t_Jtr(:,i) = t_Jtr_temp;
    else
        t_Jtr(:,i) = meshgen_t(par_pulseKeepOn);
        Jtr(:,i) = interp1(t_Jtr_temp, Jtr_temp, t_Jtr(:,i), 'spline');
    end
end

sdpsol.t_Jtr = t_Jtr;
sdpsol.Jtr = Jtr;
sdpsol.Jdk = Jdk;
sdpsol.tdwell_arr = tdwell_arr;
sdpsol.Vjump = Vjump;
sdpsol.bias_source = bias_source;
sdpsol.bias_int = bias_int;
sdpsol.probe_source = pulse_source;
sdpsol.probe_int = pulse_int;
sdpsol.probe_tmax = pulse_tmax;

end

function [Jtr, Jdk] = calcJtrJdk(sol_dwell, sol_pulseKeepOn)
    J_dwell = dfana.calcJ(sol_dwell);
    J_pulse = dfana.calcJ(sol_pulseKeepOn);
    Jdk = J_dwell.tot(end, end);
    % Use end point of dark current as baseline- could average if noisy
    Jtr = J_pulse.tot(:, end) - Jdk;
end

function par = lightRampUp(par, source, initial_int, final_int)
    switch source
        case 1
            par.g1_fun_type = 'sweep';
            % COEFF = [A_start, A_end, sweep_duration]
            par.g1_fun_arg = [initial_int, final_int, par.tmax];
        case 2
            par.g2_fun_type = 'sweep';
            % COEFF = [A_start, A_end, sweep_duration]
            par.g2_fun_arg = [initial_int, final_int, par.tmax];
    end
end

function par = lightConstant(par, source, int)
    switch source
        case 1
            par.g1_fun_type = 'constant';
            par.g1_fun_arg = int;
            if isprop(par, 'int1')
                par.int1 = int;
            end
        case 2
            par.g2_fun_type = 'constant';
            par.g2_fun_arg = int;
            if isprop(par, 'int2')
                par.int2 = int;
            end
    end
end

function sol_post = runDfLoosenTolAndReduceMaxStepFactor(sol_pre, par)
    sol_post = df(sol_pre, par);
    originalMaxStepFactor = par.MaxStepFactor;
    originalRelTol = par.RelTol;
    originalAbsTol = par.AbsTol;
    if size(sol_post.u,1) ~= par.tpoints
        % try to force smaller time steps
        par.MaxStepFactor = par.MaxStepFactor/10;
        disp([mfilename ' - the solver did not succeed, decreasing the maximum time step changing par.MaxStepFactor to ' num2str(par.MaxStepFactor)])
        sol_post = df(sol_pre, par);
    end
    if size(sol_post.u,1) ~= par.tpoints
        previousRelTol = par.RelTol;
        previousAbsTol = par.AbsTol;
        par.RelTol = min(previousRelTol*3.5, 1e-2);
        par.AbsTol = min(previousAbsTol*3.5, 1e-5);
        if par.RelTol ~= previousRelTol || par.AbsTol ~= previousAbsTol
            disp([mfilename ' - the solver did not succeed, loosening the relative tolerance from ' num2str(previousRelTol) ' to ' num2str(par.RelTol) ' and absolute tolerance from ' num2str(previousAbsTol) ' to ' num2str(par.AbsTol)])
            sol_post = df(sol_pre, par);
        end
    end
    if size(sol_post.u,1) ~= par.tpoints
        % try to force smaller time steps
        par.MaxStepFactor = par.MaxStepFactor/10;
        disp([mfilename ' - the solver did not succeed, decreasing the maximum time step changing par.MaxStepFactor to ' num2str(par.MaxStepFactor)])
        sol_post = df(sol_pre, par);
    end
    if size(sol_post.u,1) ~= par.tpoints
        previousRelTol = par.RelTol;
        previousAbsTol = par.AbsTol;
        par.RelTol = min(previousRelTol*3.5, 1e-2);
        par.AbsTol = min(previousAbsTol*3.5, 1e-5);
        if par.RelTol ~= previousRelTol || par.AbsTol ~= previousAbsTol
            disp([mfilename ' - the solver did not succeed, loosening the relative tolerance from ' num2str(previousRelTol) ' to ' num2str(par.RelTol) ' and absolute tolerance from ' num2str(previousAbsTol) ' to ' num2str(par.AbsTol)])
            sol_post = df(sol_pre, par);
        end
    end
    if size(sol_post.u,1) ~= par.tpoints
        % try to force smaller time steps
        par.MaxStepFactor = par.MaxStepFactor/10;
        disp([mfilename ' - the solver did not succeed, decreasing the maximum time step changing par.MaxStepFactor to ' num2str(par.MaxStepFactor)])
        sol_post = df(sol_pre, par);
    end
    if size(sol_post.u,1) ~= par.tpoints
        warning(['Driftfusion:' mfilename],...
            [mfilename ' - cannot make the simulation work in any way.'])
    end
    sol_post.par.MaxStepFactor = originalMaxStepFactor;
    sol_post.par.RelTol = originalRelTol;
    sol_post.par.AbsTol = originalAbsTol;
end

%------------- END OF CODE --------------