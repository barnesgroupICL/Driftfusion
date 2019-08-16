%% A script to obtain the transient ideality factor
% clear all
% close all

%% Insert the appropriate file path here
par_nid = pc('Input_files/3_layer_test.csv');

%% Get equilibrium
soleq_nid = equilibrate(par_nid);

%% Preconditioning at Vbi+0.2 V until stable- change the second argument for a different precondition
%           jumptoV(sol_ini, Vjump, tdwell, mobseti, Int, stabilise, accelerate)
sol_jump = jumptoV(soleq_nid.ion, par_nid.Vbi+0.2, 1e-3, 1, 0, 1, 1);

%% Get open circuit solutions
sol_OC = transient_nid(sol_jump, logspace(-3,2,6), 10, 1, 1e6, 200);

%% Analyse the solution
nidt_mean = transient_nid_ana(sol_OC);

% Plot initial OC EL diagram at 1 Sun
dfplot.ELx(sol_OC(5),0);
