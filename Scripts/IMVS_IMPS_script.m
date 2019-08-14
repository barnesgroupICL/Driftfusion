%% Perform IMVS measurement followed by IMPS
% clear all
% close all

par_3l = pc('Input_files/3_layer_test_symmetric.csv');

%% Get equilibrium
soleq_3l = equilibrate(par_3l);

%% Do IMPS
% do_IMVS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
sol_IMVS = do_IMVS(soleq_3l.ion, 1, 0.2, 4, 1, 200);

%% Do IMPS
sol_IMPS = do_IMPS(soleq_3l.ion, 0, 0.02, 4, 1, 200);

%% Plot the generation rate as a function of x and t
x_ihalf = getvarihalf(sol_IMVS.x);

figure(200)
surf(x_ihalf, sol_IMVS.t, sol_IMVS.g)
xlabel('Position [cm]')
ylabel('Time [s]')
zlabel('Generation rate [cm-3s-1]')