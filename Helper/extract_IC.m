function sol_ic = extract_IC(sol, requested_time)

index = find(sol.t <= requested_time);
index = index(end);

sol_ic = sol;
sol_ic.u = sol.u(index, :, :);

% Overwrite Vapp
Vappt = dfana.calcVapp(sol);
sol_ic.par.Vapp = Vappt(index);
sol_ic.t = 0;

end