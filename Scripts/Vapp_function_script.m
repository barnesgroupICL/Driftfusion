% Creates a single carrier device and then applies a 20 mV periodic
% potential for 2 cycles

par.singlecar = pc('Input_files/1_layer_single_carrier.csv');

soleq = equilibrate(par.singlecar);

% tmax is the period (seconds)
tmax = 1e-2;
Nperiods = 4;   % Number of periods
%coeff = [20e-3, Nperiods*(2*pi)/tmax,0];
coeff = [1, Nperiods*(2*pi)/tmax,0];
Vapp_func = @(coeff, t) coeff(1)*sin(coeff(2)*t + coeff(3));

% Vapp_function(sol_ini, Vapp_func, tmax, tpoints, logtime)
sol_Vapp_func = Vapp_function(soleq.ion, Vapp_func, coeff, tmax, 200, 0);

% Plot outputs
dfplot.Vappt(sol_Vapp_func)
% Current at mid-point
dfplot.Jt(sol_Vapp_func, round(par.singlecar.pcum(end)/2))
% JV plot
dfplot.JVapp(sol_Vapp_func, round(par.singlecar.pcum(end)/2))
% Energy level diagrams at t=0 and max amplitude
dfplot.ELx(sol_Vapp_func, [0, tmax/(4*Nperiods)]);

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder size regardless of whether it contains data or not.
% save('/Users/Username/Data/temp.mat')