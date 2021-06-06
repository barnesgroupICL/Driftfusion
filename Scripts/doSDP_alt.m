function sdpsol = doSDP_alt(sol_ini, tdwell_arr, Vjump, bias_source, bias_int, pulse_source, pulse_int, pulse_tmax, pulse_mobile_ions, solver_split_pulse)
%DOSDP_ALT - an alternative script for Step-Dwell-Probe simulations.
%
% Syntax:  sdpsol = doSDP_alt(sol_ini, tdwell_arr, Vjump, bias_source, bias_int, pulse_source, pulse_int, pulse_tmax, pulse_mobile_ions, solver_split_pulse)
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
%   SOLVER_SPLIT_PULSE - a boolean, whether is ok to solve the pulse
%     (probe) step breaking the time axis in smaller simulations. It can
%     help the cases in where the solver breaks when simulating this step.
%     The only drawback in activating it is that the simulation will take
%     longer.
%
% Outputs:
%   SDPSOL - a structure, including: Jdk, the current value at the end of the dwell
%     step; Jtr, the current profile during the pulse step minus the
%     Jdk current; t_Jtr, the time axis during the pulse step; and a few
%     parameters as supplied in the input arguments.
%
% Example:
%   sdpsol_spiro_01sun = doSDP_alt(soleq_spiro.ion, [logspace(-8,3,56)], 0.6, 1, 0.1, 2, 5.12, 1e-3, true, true);
%     perform 56 SDP simulations using dwell times from 10 ns to 1000 s,
%     using a voltage step of 0.6 V, illuminating from the primary light
%     source (AM15 by default) with intensity of 0.1 sun during the dwell
%     step, adding a secondary illumination (a red laser by default) with
%     an intensity of 5.12 as the pulse with a duration of 1 ms, which is
%     too long for ensuring that the ions didn't move, so that the first
%     "true" value instructs the solver to freeze the ions during this
%     step. The second "true" value allows the script to solve the pulse
%     step in smaller substeps ensuring convergence of the simulation.
%
% Other m-files required: df, dfana
% Subfunctions: none
% MAT-files required: none
%
% See also df, SDP_script_ana, export_SDP_results, doSDP.
%
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
sol_jump = runDf(sol_ini, par_jump);

par_dwell = lightConstant(par_jump, bias_source, bias_int);
% unset Vjump parameters
par_dwell.V_fun_type = 'constant';
par_dwell.V_fun_arg = Vjump;
% thaw ions
par_dwell.mobseti = 1;
par_dwell.tmesh_type = 2;
par_dwell.t0 = jump_time/10;

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

startFromSolJump = true;

for i = 1:length(tdwell_arr)
    disp(['SDP tdwell = ' num2str(tdwell_arr(i))])

    if ~startFromSolJump
        par_dwell.tmax = tdwell_arr(i) - tdwell_arr(i-1);
        sol_dwell = runDf(sol_dwell, par_dwell);
        if size(sol_dwell.u,1) ~= par_dwell.tpoints
            startFromSolJump = true;
        end
    end
    if startFromSolJump
        par_dwell.tmax = tdwell_arr(i);
        sol_dwell = runDf(sol_jump, par_dwell);
        %to save the dwell solutions (via assignin) we want to keep
        %startFromSolJump = true, so you'll have to comment out the whole
        %if/else block where it could be set to true
        %assignin('base', matlab.lang.makeValidName(['SDP_dwell_' inputname(1) '_int_' num2str(bias_int) '_tdwell_' num2str(tdwell_arr(i))]), sol_dwell)
        if size(sol_dwell.u,1) == par_dwell.tpoints
            startFromSolJump = false;
        else
            sol_dwell = runDfFourTrunks(sol_jump, par_dwell, '');
            % this block is for saving the dwell solution
%             [sol_dwell, sol1, sol2, sol3] = runDfFourTrunks(sol_jump, par_dwell, '');
%             J1 = dfana.calcJ(sol1);
%             J2 = dfana.calcJ(sol2);        
%             J3 = dfana.calcJ(sol3);
%             J4 = dfana.calcJ(sol_dwell);
%             Jtot1 = J1.tot(:,end);
%             Jtot2 = J2.tot(:,end);        
%             Jtot3 = J3.tot(:,end);
%             Jtot4 = J4.tot(:,end);
% 
%             Jtot_dwell = [Jtot1; Jtot2; Jtot3; Jtot4];
%             size(Jtot_dwell)
%             t_Jtot_dwell = [sol1.t'; (sol2.t+sol1.par.tmax)'; (sol3.t+sol1.par.tmax+sol2.par.tmax)'; (sol_dwell.t+sol1.par.tmax+sol2.par.tmax+sol3.par.tmax)'];
%             size(t_Jtot_dwell)
%             assignin('base', matlab.lang.makeValidName(['SDP_dwell_Jtot_' inputname(1) '_int_' num2str(bias_int) '_tdwell_' num2str(tdwell_arr(i))]), [t_Jtot_dwell, Jtot_dwell])
            if size(sol_dwell.u,1) == par_dwell.tpoints
                startFromSolJump = false;
            end
        end
    end
    
    disp('Pulse turn on')
    sol_pulseTurnOn = runDf(sol_dwell, par_pulseTurnOn);
    
    disp('Pulse keep on')
    if solver_split_pulse
        [~, ~, ~, ~, Jtr(:,i), Jdk(i), t_Jtr(:,i)]= runDfFourTrunks(sol_pulseTurnOn, par_pulseKeepOn, '', sol_dwell);
    else
        sol_pulseKeepOn = runDf(sol_pulseTurnOn, par_pulseKeepOn);
        [Jtr(:,i), Jdk(i)] = calcJtrJdk(sol_dwell, sol_pulseKeepOn);
        t_Jtr(:,i) = sol_pulseKeepOn.t;
        %assignin('base', matlab.lang.makeValidName(['SDP_pulseKeepOn' inputname(1) '_' num2str(tdwell_arr(i))]), sol_pulseKeepOn)
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

function sol_post = runDf(sol_pre, par)
    sol_post = df(sol_pre, par);
    originalMaxStepFactor = par.MaxStepFactor;
    originalRelTol = par.RelTol;
    originalAbsTol = par.AbsTol;
    if size(sol_post.u,1) ~= par.tpoints
        % try to force smaller time steps
        par.MaxStepFactor = par.MaxStepFactor/10;
        %warning(['Driftfusion:' mfilename],...
        disp([mfilename ' - the solver did not succeed, decreasing the maximum time step changing par.MaxStepFactor to ' num2str(par.MaxStepFactor)])
        sol_post = df(sol_pre, par);
    end
    if size(sol_post.u,1) ~= par.tpoints
        previousRelTol = par.RelTol;
        previousAbsTol = par.AbsTol;
        par.RelTol = min(previousRelTol*3.5, 1e-2);
        par.AbsTol = min(previousAbsTol*3.5, 1e-5);
        if par.RelTol ~= previousRelTol || par.AbsTol ~= previousAbsTol
            %warning(['Driftfusion:' mfilename],...
            disp([mfilename ' - the solver did not succeed, loosening the relative tolerance from ' num2str(previousRelTol) ' to ' num2str(par.RelTol) ' and absolute tolerance from ' num2str(previousAbsTol) ' to ' num2str(par.AbsTol)])
            sol_post = df(sol_pre, par);
        end
    end
    if size(sol_post.u,1) ~= par.tpoints
        % try to force smaller time steps
        par.MaxStepFactor = par.MaxStepFactor/10;
        %warning(['Driftfusion:' mfilename],...
        disp([mfilename ' - the solver did not succeed, decreasing the maximum time step changing par.MaxStepFactor to ' num2str(par.MaxStepFactor)])
        sol_post = df(sol_pre, par);
    end
    if size(sol_post.u,1) ~= par.tpoints
        previousRelTol = par.RelTol;
        previousAbsTol = par.AbsTol;
        par.RelTol = min(previousRelTol*3.5, 1e-2);
        par.AbsTol = min(previousAbsTol*3.5, 1e-5);
        if par.RelTol ~= previousRelTol || par.AbsTol ~= previousAbsTol
            %warning(['Driftfusion:' mfilename],...
            disp([mfilename ' - the solver did not succeed, loosening the relative tolerance from ' num2str(previousRelTol) ' to ' num2str(par.RelTol) ' and absolute tolerance from ' num2str(previousAbsTol) ' to ' num2str(par.AbsTol)])
            sol_post = df(sol_pre, par);
        end
    end
    if size(sol_post.u,1) ~= par.tpoints
        % try to force smaller time steps
        par.MaxStepFactor = par.MaxStepFactor/10;
        %warning(['Driftfusion:' mfilename],...
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

function [sol4, sol1, sol2, sol3, Jtr, Jdk, t_Jtr] = runDfFourTrunks(sol_pre, par, extraOutput, sol_dwell)
        par1 = par;
        par1.tmax = (par.t0^1.5*par.tmax^0.5)^0.5;
        disp(['Split solution, part ' extraOutput '1. tmax: ' num2str(par1.tmax) ' s'])
        par1.tpoints = max(floor(par.tpoints / 4), 10);
        sol1 = runDf(sol_pre, par1);
        if size(sol1.u,1) ~= par1.tpoints
            newExtraOutput = [extraOutput '1.'];
            [sol_new4, sol_new1, sol_new2, sol_new3] = runDfFourTrunks(sol_pre, par1, newExtraOutput);
            sol1 = mergeFourSolutions(sol_new1, sol_new2, sol_new3, sol_new4);
        end

        par2 = par;
        par2.tmax = (par.t0*par.tmax)^0.5 - par1.tmax;
        disp(['Split solution, part ' extraOutput '2. tmax: ' num2str(par2.tmax) ' s'])
        par2.tpoints = max(floor(par.tpoints / 4), 10);
        par2.t0 = par1.tmax / 10;
        sol2 = runDf(sol1, par2);
        if size(sol2.u,1) ~= par2.tpoints
            newExtraOutput = [extraOutput '2.'];
            [sol_new4, sol_new1, sol_new2, sol_new3] = runDfFourTrunks(sol1, par2, newExtraOutput);
            sol2 = mergeFourSolutions(sol_new1, sol_new2, sol_new3, sol_new4);
        end
        
        par3 = par;
        par3.tmax = (par.t0^0.5*par.tmax^1.5)^0.5 - par1.tmax - par2.tmax;
        disp(['Split solution, part ' extraOutput '3. tmax: ' num2str(par3.tmax) ' s'])
        par3.tpoints = max(floor(par.tpoints / 4), 10);
        par3.t0 = par2.tmax / 10;
        sol3 = runDf(sol2, par3);
        if size(sol3.u,1) ~= par3.tpoints
            newExtraOutput = [extraOutput '3.'];
            [sol_new4, sol_new1, sol_new2, sol_new3] = runDfFourTrunks(sol2, par3, newExtraOutput);
            sol3 = mergeFourSolutions(sol_new1, sol_new2, sol_new3, sol_new4);
        end
        
        par4 = par;
        par4.tmax = par.tmax - par1.tmax - par2.tmax - par3.tmax;
        disp(['Split solution, part ' extraOutput '4. tmax: ' num2str(par4.tmax) ' s'])
        par4.tpoints = max(par.tpoints - par1.tpoints - par2.tpoints - par3.tpoints, 10);
        par4.t0 = par3.tmax / 10;
        sol4 = runDf(sol3, par4);
        if size(sol4.u,1) ~= par4.tpoints
            newExtraOutput = [extraOutput '4.'];
            [sol_new4, sol_new1, sol_new2, sol_new3] = runDfFourTrunks(sol3, par4, newExtraOutput);
            sol4 = mergeFourSolutions(sol_new1, sol_new2, sol_new3, sol_new4);
        end
        
        if nargout > 4
            [Jtr_temp1, Jdk] = calcJtrJdk(sol_dwell, sol1);
            Jtr_temp2 = calcJtrJdk(sol_dwell, sol2);        
            Jtr_temp3 = calcJtrJdk(sol_dwell, sol3);
            Jtr_temp4 = calcJtrJdk(sol_dwell, sol4);

            Jtr = [Jtr_temp1; Jtr_temp2; Jtr_temp3; Jtr_temp4];
            t_Jtr = [sol1.t'; (sol2.t+par1.tmax)'; (sol3.t+par1.tmax+par2.tmax)'; (sol4.t+par1.tmax+par2.tmax+par3.tmax)'];
        end
end

function sol = mergeFourSolutions(sol1, sol2, sol3, sol4)
    sol.u = [sol1.u; sol2.u; sol3.u; sol4.u];
    sol.x = sol4.x;
    sol.par = sol4.par;
    sol.t = [sol1.t, sol2.t+sol1.t(end), sol3.t+sol1.t(end)+sol2.t(end), sol4.t+sol1.t(end)+sol2.t(end)+sol3.t(end)];
end
