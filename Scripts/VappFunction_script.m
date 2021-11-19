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
par_singlecar = pc('Input_files/1_layer_single_carrier.csv');

soleq_singlecar = equilibrate(par_singlecar);

% Extract parameters
par = soleq_singlecar.ion.par;

% Set upi time mesh
par.tmesh_type = 1;
par.t0 = 0;
par.tmax = 1e-2;
par.tpoints = 200;

%% Define the voltage function
par.V_fun_type = 'sin';         % Voltage function type
par.V_fun_arg(1) = 0;           % DC offset voltage (V)
par.V_fun_arg(2) = 20e-3;       % AC voltage amplitude (V)
par.V_fun_arg(3) = 1e3;         % Frequency (Hz)
par.V_fun_arg(4) = 0;           % Phase (Rads)

disp('Applying oscillating potential')
sol_Vapp = df(soleq_singlecar.ion, par);

% Plot outputs
dfplot.Vappt(sol_Vapp)
% Current at mid-point
dfplot.Jt(sol_Vapp, par_singlecar.dcum(end)/2);
%ylim([-2e-4, 2e-4])
% JV plot
dfplot.JVapp(sol_Vapp, par_singlecar.dcum(end)/2);
% Energy level diagrams at t=0 and max amplitude
dfplot.ELxnpxacx(sol_Vapp, 0);

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder size regardless of whether it contains data or not.
% save('/Users/Username/Data/temp.mat')