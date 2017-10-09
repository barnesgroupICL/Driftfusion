function sol_Int = changeLight(sol, newInt)
% stabilize at a new light intensity
p = sol.params;
p.pulseon = 0;

warning('off','pindrift:verifyStabilization');
if newInt ~= p.Int % verify that the solution is not already at that light intensity
    disp(['changeLight - Go to new light intensity: ' num2str(newInt)])
    p.Int = newInt; % set new light intensity
    sol_Int = pindrift(sol, p); % stabilize at new light intensity
else
    disp('changeLight - Already at requested light intensity')
    sol_Int = sol; % otherwise just answer the input solution
end

while ~verifyStabilization(sol_Int.sol, sol_Int.t, 0.01) % check in a strict way
    disp(['changeLight - Stabilizing over ' num2str(p.tmax) ' s']);
    sol_Int = pindrift(sol_Int, p);
    p.tmax = p.tmax*2; % this value doesn't get saved in sol_Int.params.tmax unless an additional pindrift is run
end
warning('on','pindrift:verifyStabilization');


%sol_Int.params.tmax = previous_tmax; % restore tmax to original value