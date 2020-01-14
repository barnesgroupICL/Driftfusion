function sol_CV = doCV(sol_ini, light_intensity, Vmax, Vmin, scan_rate, cycles, tpoints)
% Performs a cyclic voltammogram (CV) simulation
% Input arguments:
% SOL_INI = solution containing intitial conditions
% LIGHT_INTENSITY = Light intensity for bias light (Suns)
% VMAX = Maximum voltage point (V)
% VMIN = Minimum voltage point (V)
% SCAN_RATE = Scan rate (Vs-1)
% CYCLES = No. of scan cycles
% TPOINTS = No. of points in output time array
% P. Calado, 2020, Imperial College London
par = sol_ini.par;
V0 = par.Vapp;

if light_intensity > 0
    sol = lightonRs(sol_ini, light_intensity, -1, 0, 0, 10);
else
    sol = sol_ini;
end

% Calculate tmax from scan rate and voltage points
deltaV = abs(Vmax-Vmin)+abs(Vmin-Vmax)+abs(V0-Vmin);
tmax = (scan_rate*deltaV);

disp('Performing cyclic voltamagram')
sol_CV = VappFunction(sol, 'tri', [V0, Vmax, Vmin, cycles, tmax/cycles], tmax, tpoints, 0);
disp('Complete')

end