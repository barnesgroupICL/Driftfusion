% Single carrier device test
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
%% Load parameters
par = pc('Input_files/1_layer_single_carrier.csv');

%% Run to equilibrium
soleq_single_carrier = equilibrate(par);

%% Do Cyclic Voltammograms (CV) at 1, 10, and 100 Vs-1 for 4 cycles each
% doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_k0p01 = doCV(soleq_single_carrier.ion, 0, 0, 0.6, -0.6, 1e-2, 4, 401);
sol_CV_k0p1 = doCV(soleq_single_carrier.ion, 0, 0, 0.6, -0.6, 1e-1, 4, 401);
sol_CV_k1 = doCV(soleq_single_carrier.ion, 0, 0, 0.6, -0.6, 1, 4, 401);

%% Plot CV using currents from left-hand boundary (x=0)
dfplot.JVapp(sol_CV_k0p01, 0);
hold on
dfplot.JVapp(sol_CV_k0p1, 0);
hold on
dfplot.JVapp(sol_CV_k1, 0);