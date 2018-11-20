function relax = jumptoV(sol_ini, Vjump)

par = sol_ini.p;
V0 = par.Vapp;
V1 = Vjump;

jump1 = doJV(sol_ini,1,100,0,0,V0,V1,1);

par.Vapp = V1;
par.mobseti = 1;
par.tmax = 10;
par.t0 = 0;%par.tmax/1e6;
par.tmesh_type = 1;
par.tpoints = 200;

sol = pindrift(jump1.dk.f, par);

all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

% loop to check ions have reached stable config- if not accelerate ions by
% order of mag
%{
while any(all_stable) == 0
    
    % use mobseti as a multiplier- bit dangerous? But fine for now
    par.tmax = par.tmax*10;
    %par.t0 = par.tmax/1e6;
    
    sol = pindrift(jump1.dk.f, par);
    
    all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

end
%}
relax = sol;

end