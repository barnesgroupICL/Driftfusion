% Creates a single carrier device and then applies a 20 mV periodic
% potential for 2 cycles
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
par_SJ = pc('./PEM_workshop_Input_files/doped_Shottky_junction.csv');

soleq_SJ = equilibrate(par_SJ);

%%
sol_relax = jumptoV(soleq_SJ.el, 0.6, 1e-8, 0, 0, 0, 0);

%% Extract parameters
par = soleq_SJ.el.par;

% Plot outputs
dfplot.Vappt(sol_Vapp)
% Current at mid-point
dfplot.Jt(sol_Vapp, par_SJ.d_mida);
%ylim([-2e-4, 2e-4])
% JV plot
dfplot.JVapp(sol_Vapp, par_SJ.dcum(end)/2);
% Energy level diagrams at t=0 and max amplitude
dfplot.ELxnpxacx(sol_Vapp, 0);

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder size regardless of whether it contains data or not.
% save('/Users/Username/Data/temp.mat')