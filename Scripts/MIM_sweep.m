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
par_mim = pc('Input_files/mim.csv');
soleq_mim = equilibrate(par_mim);

%% Cyclic voltammogram
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV = doCV(soleq_mim.el, 0, 0, 0.4, -0.4, 1e5, 1, 401);

%% Plots
dfplot.Jt(sol_CV, par_mim.dcum(end)/2);

dfplot.JVapp(sol_CV, par_mim.dcum(end)/2);
% Energy level diagrams at t=0 and max amplitude
dfplot.ELxnpxacx(sol_CV, 0);

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder size regardless of whether it contains data or not.
% save('/Users/Username/Data/temp.mat')