function sol_ic = extract_IC(sol, requested_time)
% EXTRACT INITIAL CONDITIONS
% This function takes a single time point from the solution SOL closest to
% REQUESTED_TIME and outputs a solution SOL_IC with only that time point. This can
% then be used as the initial conditions for a new simulation with DRIFTFUSION. The time
% array SOL_IC.T is set to zero (i.e. a single time point) and VAPP is
% calculated from the SOL parameters and REQUESTED_TIME.
% Update 02/09/21 REQUESTED_TIME can now be a 2 element vector [START_TIME,
% END_TIME]
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% START CODE
if length(requested_time) == 1
    index = 0;
    index = find(sol.t <= requested_time);
    index = index(end);
elseif length(requested_time) == 2
    index = [0,0];
    index_temp1 = find(sol.t <= requested_time(1));
    index(1) = index_temp1(end);
    index_temp2 = find(sol.t <= requested_time(2));
    index(2) = index_temp2(end);
end

sol_ic = sol;
sol_ic.u = sol.u(index, :, :);

% Overwrite Vapp
Vappt = dfana.calcVapp(sol);
sol_ic.par.Vapp = Vappt(index);
% Overwrite time array
sol_ic.t = sol.t(index);

end