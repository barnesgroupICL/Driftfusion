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

% Run to equilibrium
soleq = equilibrate(par);

% Do JV scn
JVsol = doJV(soleq.el, 100e-3, 201, 1, 0, 0, 1, 1);

% plot JV scan
dfplot.JV(JVsol, 1);
set(gca,'YScale','log')