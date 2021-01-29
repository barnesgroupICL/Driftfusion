% Perform IMVS measurement followed by IMPS
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
% clear all
% close all

%% Build the parameters object
par_3l = pc('Input_files/3_layer_test_symmetric.csv');

%% Get equilibrium solutions
soleq_3l = equilibrate(par_3l);

%% Do IMPS
% do_IMVS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
sol_IMVS = doIMVS(soleq_3l.ion, 1, 0.2, 4, 1, 401);

%% Do IMPS
sol_IMPS = doIMPS(soleq_3l.ion, 0, 0.02, 4, 1, 401);

%% Plot the Voc as a function of t
dfplot.Voct(sol_IMVS);

%% Plot the J as a function of t
dfplot.Jt(sol_IMPS,0)

%% Plot the generation rate as a function of x and t
dfplot.gxt(sol_IMVS)