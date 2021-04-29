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
% Initialise the system
initialise_df

%% Create a parameters object for Spiro/MAPI/TiO2 by including a filepath to the 
% appropriate .csv as the arugment to the parameters class PC
par = pc('Input_files/3_layer_test_vary.csv');

%% Ideal transport, no rec
par_ideal = par;
par_ideal.sn(2) = 1e-20;
par_ideal.sp(2) = 1e-20;
par_ideal.sn(4) = 1e-20;
par_ideal.sp(4) = 1e-20;
par_ideal.mue(2) = 1e6;
par_ideal.muh(2) = 1e6;
par_ideal.mue(4) = 1e6;
par_ideal.muh(4) = 1e6;
par_ideal = refresh_device(par_ideal);

soleq_ideal = equilibrate(par_ideal);

%sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_ideal = doCV(soleq_ideal.el, 1, 0, 1, 0, 100e-3, 1, 281);

%% Poor transport, no rec
par_trans = par;
par_trans.sn(2) = 1e-20;
par_trans.sp(2) = 1e-20;
par_trans.sn(4) = 1e-20;
par_trans.sp(4) = 1e-20;
par_trans.mue(2) = 1e6;
par_trans.muh(2) = 1e6;
par_trans.mue(4) = 1e6;
par_trans.muh(4) = 1e6;
par_trans = refresh_device(par_trans);

soleq_trans = equilibrate(par_trans);

%sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_trans = doCV(soleq_trans.el, 1, 0, 1, 0, 100e-3, 1, 281);