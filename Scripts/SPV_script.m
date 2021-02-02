% A tutorial script showing how to perform SPV simulations with DRIFTFUSION
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
%
%% Start code
%% Initialise Driftfusion - this puts the folders in the MATLAB file path
% Driftfusion only needs to be initialised once at the start of each
% session
initialise_df

%% Create a parameters object
% A parameters object is created by calling the parameters class PC. In
% this instance I am calling the object PAR_TIO2_2L.
% To input a different parameter set duplicate SPV_TiO2_MAPI_bilayer.csv,
% open it in Excel and then edit the parameters. Save it as a different
% name and then use that as the file in the path of the input argument to
% PC (the parameters class):
par_tio2_2l = pc('input_files/tio2_mapi_bilayer.csv');

%% Get equilibrium solutions
% You will find that with the default parameters, contained in the above
% .csv file, the hole Fermi level is not well calculated at equilibrium. This is
% because the hole densities are very small and in the current version
% of Driftfusion we solve for carrier densities as opposed to Fermi levels. 
% Since the hole densities are effectively zero in the dark it's safe to 
% ignore them since they do not contribute to currents or the electrostatic 
% potential significantly.
soleq_tio2_2l = equilibrate(par_tio2_2l);

%% Do the SPV procedure
% SOLEQ_TIO2_2L contains 2 solutions, with and without ions. Here we will use
% the solution SOLEQ_TIO2_2L.ION as the input to DOSPV which can be found in
% the PROTOCOLS folder. The code is commented with the information you need
% regarding the input arguments. The code splits the solution into 2 parts-
% the illuminated step and the dark step, SPVSOL.ILL and SPVSOL.DK.
% Input arguments: dospv(sol_ini, Int, mobseti, tpoints, tmax, Rs, stabilise)
spvsol = doSPV(soleq_tio2_2l.ion, 0.2, 1, 200, 10, 1e6, 0);

%% Analyse the solution
% SPVANA (in the Analysis folder) takes the two solutions contained in SPVSOL, extracts the charge
% density, electric field, and potential at the right hand boundary as a 
% function of time the plots them. The outputs can be store in the output 
% structure SPVDAT (although this can be called anything).
spvdat = spvana(spvsol);

%% Save your workspace
% You can rename this and move it later
%save('SPV_temp.mat')

%% Plotting inidividual solutions
% Various plots are available in the plotting class DFPLOT (You can find this in the 
% CORE folder). Some of these may not work but the main ones should.  By way of example 
% if you wanted to plot the final Energy level diagram, carrier densities and 
% ionic charge densities for the illumated step of the SPV then type:
dfplot.ELxnpxacx(spvsol.ill);

% For plots that are a function of position you can use a second argument 
% to overlay the plots for different times. In this example I plot the 
% electrostatic potential as a function of position for 5 different times
dfplot.Vx(spvsol.ill, [0,1,2,4,8]);
