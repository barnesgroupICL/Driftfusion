function sol_OC = transient_nid(sol_ini, int_arr, stab_time, mobseti, Rs, pnts)
% Performs a transient ideality factor simulation. 
% See Calado 2019, Physical Review Applied for details of the measurement protocol.
% Input arguments:
% SOL_INI = initial conditions
% INT_ARR = An array of light intensities (Suns)
% STAB_TIME = Stabilisation time - length of the transient. Setting to -1
% to cycle to stable solution from initial time of 0.1
% MOBSETI = Ion mobility switch
% RS = Series resistance - recommended to use Rs = 1e6 for approx open
% circuit
% PNTS = Number of time points in the solution
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%% Start code
for i = 1:length(int_arr)
    sol_OC(i) = lightonRs(sol_ini, int_arr(i), stab_time, mobseti, Rs, pnts);
end

end