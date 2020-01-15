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
% P. Calado, 2019, Imperial College London

for i = 1:length(int_arr)
    sol_OC(i) = lightonRs(sol_ini, int_arr(i), stab_time, mobseti, Rs, pnts);
end

end