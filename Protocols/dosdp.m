function sdpsol = dosdp(sol_ini, tdwell_arr, Vjump, pulselen, pulseint, scalefactor)

% tarr =  a time array containing the dwell times
% Vjump = the jump to voltage. Vpre is defined by sol_ini which should be
% a solution at steady-state
% pulselen = length of pulse
% pulseint = pulse intensity
% Jtr_time = time at which is the J value is taken
% scalefactor = accelerates ions in case of bad convergence. Set to

% Fast initial scan
par = sol_ini.par;
t = sol_ini.t;
V0 = par.Vapp;
V1 = Vjump;

jump1 = doJV(sol_ini,1,100,0,0,V0,V1,1);

par.t0 = 0;
par.tmesh_type = 1;
par.tpoints = 100;
par.Vapp = V1;
par.mobseti = 1;
Jtr = zeros(length(t), length(tdwell_arr));

% 1 ns stabilisation
par.int1 = 0;
par.Vapp = Vjump;
par.JV = 0;

sol_stab = df(jump1.dk.f, par);

for i = 1:length(tdwell_arr)
    msg = ['SDP scan, ', num2str(i), ' tdwell = ', num2str(tdwell_arr(i))];
    disp(msg)
    
    par.int1 = 0;
    par.tmax = tdwell_arr(i)./scalefactor;
    %accelerate ions
    par.mobseti = scalefactor;
    
    sol_dk = df(sol_stab, par);
    
    % Use end point of dark current as baseline- could average if noisy
    sdpsol.Jdk = sol_dk.Jtotr(end, end);
    
    par.tmax = pulselen;
    par.int1 = pulseint;
    % switch off ion during pulse
    par.mobseti = 1;
    
    sol_pulse = df(sol_dk, par);
    sdpsol.t_Jtr = sol_pulse.t;
    
    % Jtr has baseline already removed
    sdpsol.Jtr(:,i) = sol_pulse.Jtotr(:, end) - sdpsol.Jdk;
    
end

sdpsol.tdwell_arr = tdwell_arr*scalefactor;
sdpsol.Vjump = Vjump;

sdpanal(sdpsol, 4e-6)

end





