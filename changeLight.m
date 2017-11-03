function sol_Int = changeLight(sol, newInt, tmax)
% stabilize at a new light intensity
p = sol.params;
p.pulseon = 0;
p.tmesh_type = 2;
p.tmax = tmax;
p.t0 = p.tmax / 1e11;
p.tpoints = 100;

warning('off','pindrift:verifyStabilization'); % warnings about stability are not needed here
% even if the solution is already at the requested light intensity,
% stabilize it again
disp(['changeLight - Go from light intensity ' num2str(p.Int) ' to ' num2str(newInt)])
p.Int = newInt; % set new light intensity
sol_Int = pindrift(sol, p); % stabilize at new light intensity

while ~verifyStabilization(sol_Int.sol, sol_Int.t, 0.000001) % check stability in a strict way
    disp(['changeLight - Stabilizing over ' num2str(p.tmax) ' s']);
    sol_Int = pindrift(sol_Int, p);
    p.tmax = p.tmax * 10; % this new value doesn't get saved in sol_Int.params.tmax unless an additional pindrift is run
end
warning('on','pindrift:verifyStabilization');

% just repeat the last one, for sake of paranoia yeee
sol_Int = pindrift(sol_Int, sol_Int.params);