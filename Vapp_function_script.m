% Creates a single carrier device and then applies a 50 mV periodic potential

par.singlecar = pc('input_files/1 layer single carrier.csv');

soleq = equilibrate(par.singlecar);

% tmax is the period (seconds)
tmax = 10;
coeff = [50e-3,(2*pi)/tmax,0];
Vapp_func = @(coeff, t) coeff(1)*sin(coeff(2)*t + coeff(3));

% Vapp_function(sol_ini, Vapp_func, tmax, tpoints, logtime)
sol_Vapp_func = Vapp_function(soleq.ion, Vapp_func, coeff, tmax, 200, 0);

% Plot outputs
dfplot.Vappt(sol_Vapp_func)
% Current at mid-point
dfplot.Jt(sol_Vapp_func, round(par.singular.pcum(end)/2))
