% Single layer test script
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
% Load parameters
par = pc('Input_files/1_layer_test.csv');

%% Equilibrium oslutions
soleq = equilibrate(par);
% Find equilibrium solution with 100 Ohms series resistance in external
% circuit
sol_Rs = lightonRs(soleq.el, 0, -1e-3, 0, 100, 20);

%% Do JV scans
JVsol_Rs0 = doJV(soleq.el, 100e-3, 201, 1, 0, 0, 1, 1);
JVsol_Rs100 = doJV(sol_Rs, 100e-3, 201, 1, 0, 0, 1, 1);

%% plot JV scan
dfplot.JV(JVsol_Rs0, 1);
hold on
dfplot.JV(JVsol_Rs100, 1);
hold off
set(gca,'YScale','log')