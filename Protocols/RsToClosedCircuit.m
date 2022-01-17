function sol_closed = RsToClosedCircuit(sol_Rs)
% Takes a solution with high Rs boundary condition, switches off the series
% resistance and applies Vapp equal to whatever the series resistance
% voltage was. This results in a more stable solution than that obtained
% using findDirectVoc for example.
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
par = sol_Rs.par;

% Get QFL splitting at electrode of SOL_RS
QFLsplit_t = dfana.calcDeltaQFL(sol_Rs);
QFLsplit = QFLsplit_t(end);

par.Rs = 0;
par.Vapp = QFLsplit;
% Characteristic diffusion time
t_diff = (par.dcum0(end)^2)/(2*par.kB*par.T*min(min(par.mu_n), min(par.mu_p)));
par.tmax = t_diff;
par.t0 = par.tmax/1e6;
par.tmesh_type = 2;

sol_closed = df(sol_Rs, par);

end