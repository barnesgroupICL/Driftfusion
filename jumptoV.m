function sol_relax = jumptoV(sol_ini, Vjump, tdwell)
% A function to simulate a jump-to-Voltage measurement as used for SDP
%% Input arugments
% SOL_INI - input solution- could be cell at equilibrium or stablised at an
% applied voltage
% VJUMP - the voltage to be jumped to
% TDWELL - the relaxation time after jumping

%% OUTPUTS
% SOL_RELAX - the output relaxation solution

par = sol_ini.par;
V0 = par.Vapp;
V1 = Vjump;

jump1 = doJV(sol_ini,1,100,0,0,V0,V1,1);

par.Vapp = V1;
par.mobseti = 1;
par.tmax = tdwell;
par.t0 = 0;%par.tmax/1e6;
par.tmesh_type = 1;
par.tpoints = 200;

sol = df(jump1.dk.f, par);

all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

% loop to check ions have reached stable config- if not accelerate ions by
% order of mag
%{
while any(all_stable) == 0
    
    % use mobseti as a multiplier- bit dangerous? But fine for now
    par.tmax = par.tmax*10;
    %par.t0 = par.tmax/1e6;
    
    sol = df(jump1.dk.f, par);
    
    all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

end
%}
sol_relax = sol;

end