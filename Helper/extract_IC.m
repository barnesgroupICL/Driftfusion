function sol_ic = extract_IC(sol, requested_time)
% EXTRACT INITIAL CONDITIONS
% This function takes a single time point from the solution SOL closest to
% REQUESTED_TIME and outputs a solution SOL_IC with only that time point. This can
% then be used as the initial conditions for a new simulation with DRIFTFUSION. The time
% array SOL_IC.T is set to zero (i.e. a single time point) and VAPP is
% calculated from the SOL parameters and REQUESTED_TIME.
% 
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% START CODE
index = find(sol.t <= requested_time);
index = index(end);

sol_ic = sol;
sol_ic.u = sol.u(index, :, :);

% Overwrite Vapp
Vappt = dfana.calcVapp(sol);
sol_ic.par.Vapp = Vappt(index);
% Overwrite time array
sol_ic.t = 0;

end