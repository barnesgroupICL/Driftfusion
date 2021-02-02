% This scripts demonstrates the use of the REFRESH_DEVICE function to
% rebuild the device strutcures PAR.DEV, PAR.DEV_IHALF, and the grids PAR.X and
% PAR.X_IHALF.
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
%% Initialise Driftfusion
initialise_df

%% Build a base parameter set object
par_tio2_base = pc('Input_files/spiro_mapi_tio2.csv');

%% copy to a new object where we will chaneg the recombination parameters
par_tio2_low_rec = par_tio2_base;

%% Change the interfacial recombination parameters for interface 2 (layer no. 4)
par_tio2_low_rec.sn(4) = 1;
par_tio2_low_rec.sp(4) = 1;

%% Refresh the device- rebuilds the device structures par.dev and par.dev_ihalf etc.
par_tio2_low_rec = refresh_device(par_tio2_low_rec);

%% Get equilibrium solutions
soleq_tio2_base = equilibrate(par_tio2_base);
soleq_tio2_low_rec = equilibrate(par_tio2_low_rec);

%% DoJVs 0 - 1.2 V
JV_tio2_base = doJV(soleq_tio2_base.ion, 200e-3, 200, 1, 1, 0, 1.2, 3);
JV_tio2_low_rec = doJV(soleq_tio2_low_rec.ion, 200e-3, 200, 1, 1, 0, 1.2, 3);

%% Plot JVs
dfplot.JV(JV_tio2_base, 3)
figure(4)
hold on
dfplot.JV(JV_tio2_low_rec, 3)
ylim([-25e-3, 10e-3])
legend('Base, dk F', 'Base, dk R', 'Base, 1sun F', 'Base, 1sun R',...
    'Low rec, dk F', 'Low rec, dk R', 'Low rec, 1sun F', 'Low rec, 1sun R')
hold off
