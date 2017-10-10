function sol_Int = changeLight(sol, newInt)
% stabilize at a new light intensity
p = sol.params;
p.pulseon = 0;

warning('off','pindrift:verifyStabilization'); % warnings about stability are not needed here
% even if the solution is already at the requested light intensity,
% stabilize it again
disp(['changeLight - Go from light intensity ' num2str(p.Int) ' to ' num2str(newInt)])
p.Int = newInt; % set new light intensity
sol_Int = pindrift(sol, p); % stabilize at new light intensity

while ~verifyStabilization(sol_Int.sol, sol_Int.t, 0.01) % check in a strict way
    disp(['changeLight - Stabilizing over ' num2str(p.tmax) ' s']);
    sol_Int = pindrift(sol_Int, p);
    p.tmax = p.tmax*2; % this value doesn't get saved in sol_Int.params.tmax unless an additional pindrift is run
end
warning('on','pindrift:verifyStabilization');


%sol_Int.params.tmax = previous_tmax; % restore tmax to original value
