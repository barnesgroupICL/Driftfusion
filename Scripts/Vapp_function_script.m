% Creates a single carrier device and then applies a 20 mV periodic
% potential for 2 cycles
% clear all
% close all

par_singlecar = pc('Input_files/1_layer_single_carrier.csv');

soleq_singlecar = equilibrate(par_singlecar);

% Extract parameters
par = soleq_singlecar.ion.par;

% Set upi time mesh
par.tmesh_type = 1;
par.t0 = 0;
par.tmax = 1e-2;
par.tpoints = 200;

par.V_fun_type = 'sin';
par.V_fun_arg(1) = 0;
par.V_fun_arg(2) = 20e-3;
par.V_fun_arg(3) = 1e3;
par.V_fun_arg(4) = 0;

disp('Applying oscillating potential')
sol_Vapp = df(soleq_singlecar.ion, par);

% Plot outputs
dfplot.Vappt(sol_Vapp)
% Current at mid-point
dfplot.Jt(sol_Vapp, round(par_singlecar.pcum(end)/2))
%ylim([-2e-4, 2e-4])
% JV plot
dfplot.JVapp(sol_Vapp, round(par_singlecar.pcum(end)/2))
% Energy level diagrams at t=0 and max amplitude
dfplot.ELx(sol_Vapp, 0);

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder size regardless of whether it contains data or not.
% save('/Users/Username/Data/temp.mat')