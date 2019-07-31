function sol_OC = transient_nid(sol_ini, int_arr, stab_time, mobseti, Rs, pnts)

% SOL_INI = initial conditions
% INT_ARR = An array of light intensities (Suns)
% STAB_TIME = Stabilisation time - length of the transient. Setting to -1
% to cycle to stable solution from initial time of 0.1
% MOBSETI = Ion mobility switch
% RS = Series resistance - recommended to use Rs = 1e6 for approx open
% circuit
% PNTS = Number of time points in the solution

for i = 1:length(int_arr)
    sol_OC(i) = lighton_Rs(sol_ini, int_arr(i), stab_time, mobseti, Rs, pnts);
end

end