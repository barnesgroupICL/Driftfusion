 function sol_Int = changeLight(sol, newInt, tmax)
% stabilize at a new light intensity

p = sol.params;
p.pulseon = 0;
p.tmesh_type = 2;
p.t0 = 1e-10;
p.tpoints = 30;
sol_Int = sol;

% set an initial time for stabilization tmax
if tmax
    p.tmax = tmax;
else
    if p.mui
        p.tmax = 2^(-log10(p.mui)) / 10 + 2^(-log10(p.mue_i));
    else
        p.tmax = 2^(-log10(p.mue_i));
    end
end

% warnings about stability are not needed in the first steps
warning('off','pindrift:verifyStabilization'); 

% change light intensity in small steps
steps = 1 + ceil(abs(log10(newInt/p.Int)));

% if the step is just one, then newInt is the output of logspace function
Int_array = logspace(log10(p.Int), log10(newInt), steps); 

% change light in steps, not needed to reach a good stabilized solution in
% each step, so stabilization is not verified here
% skip first value in the array as is the initial Int
for i = 2:length(Int_array)
    disp([mfilename ' - Go from light intensity ' num2str(p.Int) ' to ' num2str(Int_array(i)) ' over ' num2str(p.tmax) ' s'])
    p.Int = Int_array(i); % set new light intensity
    sol_Int = pindrift(sol_Int, p);
end

while ~verifyStabilization(sol_Int.sol, sol_Int.t, 1e-8) % check stability in a strict way
    disp([mfilename ' - Stabilizing over ' num2str(p.tmax) ' s']);
    sol_Int = pindrift(sol_Int, p);
    p.tmax = p.tmax * 10; % this change gets saved in sol_Int only if another cycle is performed
end

% warnings about stability can be useful again
warning('on','pindrift:verifyStabilization');

% just repeat the last one, for sake of paranoia yeee
disp([mfilename ' - Stabilizing over ' num2str(p.tmax) ' s']);
sol_Int = pindrift(sol_Int, p);
