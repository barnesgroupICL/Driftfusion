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
par = pc('Input_files/Schottky_junction.csv');
par.MaxStepFactor = 0.1; % Reduce maximum time step

%% Run to equilibrium
soleq = equilibrate(par);

dfplot.ELx(soleq.el)

%% Electronic carriers only
%% Initial ramp
deltaV_large = 0.5;
t_ramp = 1e-12;
t_dwell = 1e-6;
[sol_ramp_large, sol_dwell_large, sol_concat_large] = ramped_step(soleq.el, deltaV_large, t_ramp, t_dwell);

%% Plots
dfplot.Jt(sol_concat_large, 0);
set(gca, 'XScale','log');
dfplot.Vappt(sol_ramp_large);

%% Small perturbation
deltaV_small = 0.05;
t_ramp = 1e-12;
t_dwell = 1e-6;
[sol_ramp_small, sol_dwell_small, sol_concat_small] = ramped_step(sol_dwell_large, deltaV_small, t_ramp, t_dwell);

%% Plots
dfplot.Jt(sol_concat_small, 0);
set(gca, 'XScale','log');
dfplot.Vappt(sol_ramp_small);

%% With ions
%% Initial ramp
deltaV_large = 0.5;
t_ramp = 1e-12;
t_dwell = 100;
[sol_ramp_large_ion, sol_dwell_large_ion, sol_concat_large_ion] = ramped_step(soleq.ion, deltaV_large, t_ramp, t_dwell);

%% Plots
dfplot.Jt(sol_concat_large_ion, 0);
set(gca, 'XScale','log');
dfplot.Vappt(sol_ramp_large_ion);

%% Small perturbation
deltaV_small = 0.05;
t_ramp = 1e-12;
t_dwell = 10;
[sol_ramp_small_ion, sol_dwell_small_ion, sol_concat_small_ion] = ramped_step(soleq.ion, deltaV_small, t_ramp, t_dwell);

%% Plots
dfplot.Jt(sol_concat_small_ion, 0);
set(gca, 'XScale','log');
dfplot.Vappt(sol_ramp_small_ion);
