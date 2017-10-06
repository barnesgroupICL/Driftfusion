function sol_Int = changeLight(sol, newInt)
% stabilize at a new light intensity
p = sol.params;

if newInt ~= p.Int % verify that the solution is not already at that light intensity
    disp(['Go to new light intensity' newInt])
    p.Int = newInt; % set new light intensity
    sol_Int = pindrift(sol, p); % stabilize at new light intensity
else
    sol_Int = sol; % otherwise just answer the input solution
end
