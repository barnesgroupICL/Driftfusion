% Single layer test
clear all
close all
% Load parameters
par = pc('Input_files/1_layer_single_carrier.csv');

% Run to equilibrium
soleq = equilibrate(par);

% Do JV scn
JVsol = doJV(soleq.el, 100e-3, 201, 1, 0, 0, 1, 1);

% plot JV scan
dfplot.JV(JVsol, 1);
