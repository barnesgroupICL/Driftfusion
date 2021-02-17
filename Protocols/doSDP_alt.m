function sdpsol = doSDP_alt(sol_ini, tdwell_arr, Vjump, bias_source, bias_int, pulse_source, pulse_int, pulse_tmax, pulse_freeze_ions, solver_split_pulse)
% TARR =  a time array containing the dwell times
% VJUMP = the jump to voltage. Vpre is defined by sol_ini which should be
% a solution at steady-state
% PULSE_TMAX = capture length for transient
% PULSE_INT = pulse intensity
% Jtr_time = time at which is the J value is taken
% scalefactor = accelerates ions in case of bad convergence. Set to
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Start code
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
par_pulseKeepOn.mobseti = pulse_freeze_ions;
par_pulseKeepOn.tmesh_type = 2;
par_pulseKeepOn.t0 = jump_time;
par_pulseKeepOn.tmax = pulse_tmax;
par_pulseKeepOn.tpoints = 100;
%par_pulseKeepOn.RelTol = 1e-4;
%par_pulseKeepOn.AbsTol = 1e-7;

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
        %assignin('base', matlab.lang.makeValidName(['SDP_dwell_' inputname(1) '_int_' num2str(bias_int) '_' num2str(tdwell_arr(i))]), sol_dwell)
        if size(sol_dwell.u,1) == par_dwell.tpoints
            startFromSolJump = false;
        else
            sol_dwell = runDfFourTrunks(sol_jump, par_dwell);
            if size(sol_dwell.u,1) == par_dwell.tpoints
                startFromSolJump = false;
            end
        end
    end
    
    disp('Pulse turn on')
    sol_pulseTurnOn = runDf(sol_dwell, par_pulseTurnOn);
    
    disp('Pulse keep on')
    if solver_split_pulse
        [~, Jtr(:,i), Jdk(i), t_Jtr(:,i)]= runDfFourTrunks(sol_pulseTurnOn, par_pulseKeepOn, sol_dwell);
    else
        sol_pulseKeepOn = runDf(sol_pulseTurnOn, par_pulseKeepOn);
        [Jtr(:,i), Jdk(i)] = calcJtrJdk(sol_dwell, sol_pulseKeepOn);
        t_Jtr(:,i) = sol_pulseKeepOn.t;
        %assignin('base', matlab.lang.makeValidName(['SDP_' inputname(1) '_' num2str(tdwell_arr(i))]), sol_pulseKeepOn)
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
            warning(['Driftfusion:' mfilename],...
                [mfilename ' - the solver did not succeed, decreasing the maximum time step changing par.MaxStepFactor to ' num2str(par.MaxStepFactor)])
        sol_post = df(sol_pre, par);
    end
    if size(sol_post.u,1) ~= par.tpoints
        par.RelTol = min(par.RelTol*3.5, 1e-2);
        par.AbsTol = min(par.AbsTol*3.5, 1e-5);
        warning(['Driftfusion:' mfilename],...
            [mfilename ' - the solver did not succeed, loosening the relative and absolute tolerance.'])
        sol_post = df(sol_pre, par);
    end
    if size(sol_post.u,1) ~= par.tpoints
        % try to force smaller time steps
        par.MaxStepFactor = par.MaxStepFactor/10;
            warning(['Driftfusion:' mfilename],...
                [mfilename ' - the solver did not succeed, decreasing the maximum time step changing par.MaxStepFactor to ' num2str(par.MaxStepFactor)])
        sol_post = df(sol_pre, par);
    end
    if size(sol_post.u,1) ~= par.tpoints
        par.RelTol = min(par.RelTol*3.5, 1e-2);
        par.AbsTol = min(par.AbsTol*3.5, 1e-5);
        warning(['Driftfusion:' mfilename],...
            [mfilename ' - the solver did not succeed, loosening the relative and absolute tolerance.'])
        sol_post = df(sol_pre, par);
    end
    if size(sol_post.u,1) ~= par.tpoints
        % try to force smaller time steps
        par.MaxStepFactor = par.MaxStepFactor/10;
            warning(['Driftfusion:' mfilename],...
                [mfilename ' - the solver did not succeed, decreasing the maximum time step changing par.MaxStepFactor to ' num2str(par.MaxStepFactor)])
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

function [sol4, Jtr, Jdk, t_Jtr] = runDfFourTrunks(sol_pre, par, sol_dwell)
        par1 = par;
        par1.tmax = (par.t0^1.5*par.tmax^0.5)^0.5;
        disp(['Split solution, first part. tmax: ' num2str(par1.tmax) ' s'])
        par1.tpoints = floor(par.tpoints / 4);
        sol1 = runDf(sol_pre, par1);

        par2 = par;
        par2.tmax = (par.t0*par.tmax)^0.5 - par1.tmax;
        disp(['Split solution, second part. tmax: ' num2str(par2.tmax) ' s'])
        par2.tpoints = floor(par.tpoints / 4);
        par2.t0 = par1.tmax / 10;
        sol2 = runDf(sol1, par2);
        
        par3 = par;
        par3.tmax = (par.t0^0.5*par.tmax^1.5)^0.5 - par1.tmax - par2.tmax;
        disp(['Split solution, third part. tmax: ' num2str(par3.tmax) ' s'])
        par3.tpoints = floor(par.tpoints / 4);
        par3.t0 = par2.tmax / 10;
        sol3 = runDf(sol2, par3);
        
        par4 = par;
        par4.tmax = par.tmax - par1.tmax - par2.tmax - par3.tmax;
        disp(['Split solution, fourth part. tmax: ' num2str(par4.tmax) ' s'])
        par4.tpoints = par.tpoints - par1.tpoints - par2.tpoints - par3.tpoints;
        par4.t0 = par3.tmax / 10;
        sol4 = runDf(sol3, par4);
        
        if nargout > 1
            [Jtr_temp1, Jdk] = calcJtrJdk(sol_dwell, sol1);
            Jtr_temp2 = calcJtrJdk(sol_dwell, sol2);        
            Jtr_temp3 = calcJtrJdk(sol_dwell, sol3);
            Jtr_temp4 = calcJtrJdk(sol_dwell, sol4);

            Jtr = [Jtr_temp1; Jtr_temp2; Jtr_temp3; Jtr_temp4];
            t_Jtr = [sol1.t'; (sol2.t+par1.tmax)'; (sol3.t+par1.tmax+par2.tmax)'; (sol4.t+par1.tmax+par2.tmax+par3.tmax)'];
        end
end