%% A script to obtain the transient ideality factor

%% Insert the appropriate file path here
par = pc('input_files/3_layer_test.csv');

%% Get equilibrium
soleq = equilibrate(par);

%% Preconditioning at 1.2 V until stable
sol_relax = jumptoV(soleq.ion, par.Vbi+0.2, 1e-3, 0, 1);

%% Get open circuit solutions
sol_OC = transient_nid(sol_relax, logspace(-4,2,7), 10, 1, 1e6, 200);

%% Analyse the solution
nidt_mean = transient_nid_ana(sol_OC);