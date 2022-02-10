%% LICENSE
% Copyright (C) 2022  Philip Calado, Ilario Gelmetti, Lucy Hart and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%% Read in parameters for 3-layer perovskite solar cell
par = pc('Input_files/spiro_mapi_tio2.csv');

%% Obtain equilibrium solution
sol_eq = equilibrate(par);

%% Run current-voltage scans for different applied biases
sol_CV1 = doCV(sol_eq.el, 1, 0.7, 1.0, 0, 100e-3, 3, 241);
sol_CV2 = doCV(sol_eq.ion, 1, 0, 1.0, 0, 1e-3, 1, 241);
sol_CV3 = doCV(sol_eq.ion, 1, 0, 1.0, 0, 500e-3, 1, 241);
sol_CV4 = doCV(sol_eq.ion, 1, 0, 0.9, 0, 1000e-3, 1, 241);

%% Plot outputs
dfplot.JtotVapp(sol_CV1, 0);
hold on
dfplot.JtotVapp(sol_CV2, 0);
hold on
dfplot.JtotVapp(sol_CV3, 0);
hold on
dfplot.JtotVapp(sol_CV4, 0);
hold off
ylim([-30e-3, 10e-3])

%% Get statistics
stats1 = CVstats(sol_CV1)
stats2 = CVstats(sol_CV2)
stats3 = CVstats(sol_CV3)
stats4 = CVstats(sol_CV4)