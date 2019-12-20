function sdpsol = doSDP(sol_ini, tdwell_arr, Vjump, pulse_int, pulse_tmax, duty, tpoints, scalefactor)
% tarr =  a time array containing the dwell times
% Vjump = the jump to voltage. Vpre is defined by sol_ini which should be
% a solution at steady-state
% pulse_tmax = capture length for transient
% pulse_int = pulse intensity
% Jtr_time = time at which is the J value is taken
% scalefactor = accelerates ions in case of bad convergence. Set to
disp('Starting SDP protocol')

% Get parameters
par = sol_ini.par;

% Voltage jump
par.t0 = 0;
par.tmax = 1e-3;
par.tmesh_type = 1;
par.tpoints = 10;
par.mobseti = 0;

par.V_fun_type = 'sweep';
par.V_fun_arg(1) = par.Vapp;        % Start at input solution Vapp
par.V_fun_arg(2) = Vjump;
par.V_fun_arg(3) = par.tmax;

disp('Initial jump')
sol_jump = df(sol_ini, par);
disp('Complete')

% Dwell phase parameters
par.t0 = 0;
par.tmesh_type = 1;
par.V_fun_type = 'constant';
par.V_fun_arg = Vjump;
par.mobseti = 1;
par.K_cation = scalefactor;
par.K_anion = scalefactor;

% Preallocate memory
sdpsol.Jtr = zeros(tpoints, length(tdwell_arr));

for i = 1:length(tdwell_arr)
    msg = ['SDP scan iteration ', num2str(i), ', tdwell = ', num2str(tdwell_arr(i))];
    disp(msg)
    
    par.tmax = tdwell_arr(i);
    
    disp('Obtaining dwell solution')
    sol_dwell = df(sol_jump, par);
    disp('Complete')
    
    % Obtain basline current
    [~,J_dwell,~] = dfana.calcJ(sol_dwell);
    sdpsol.Jdk(i) = J_dwell.tot(end, end);
    
    % Use end point of dark current as baseline- could average if noisy
    sol_pulse = doLightPulse(sol_dwell, pulse_int, pulse_tmax, tpoints, duty, 0, 1);
    
    [~,J_pulse,~] = dfana.calcJ(sol_pulse);
    sdpsol.Jtr(:,i) = J_pulse.tot(:, end) - sdpsol.Jdk(i);
end

sdpsol.t_Jtr = sol_pulse.t;
sdpsol.tdwell_arr = tdwell_arr*scalefactor;
sdpsol.Vjump = Vjump;

end





