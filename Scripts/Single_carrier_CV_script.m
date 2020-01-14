% Single carrier device test

% Load parameters
par = pc('Input_files/1_layer_single_carrier.csv');

% Run to equilibrium
soleq_single_carrier = equilibrate(par);

% Do Cyclic Voltammograms (CV) at 1, 10, and 100 Vs-1 for 4 cycles each
% doCV(sol_ini, light_intensity, Vmax, Vmin, scan_rate, cycles, tpoints)
CVsol_k1 = doCV(soleq_single_carrier.ion, 0, 0.6, -0.6, 1, 4, 401);
CVsol_k10 = doCV(soleq_single_carrier.ion, 0, 0.6, -0.6, 10, 4, 401);
CVsol_k100 = doCV(soleq_single_carrier.ion, 0, 0.6, -0.6, 100, 4, 401);

% plot CV using currents from left-hand boundary (x=0)
dfplot.JVapp(CVsol_k1, 0);
hold on
dfplot.JVapp(CVsol_k10, 0);
hold on
dfplot.JVapp(CVsol_k100, 0);