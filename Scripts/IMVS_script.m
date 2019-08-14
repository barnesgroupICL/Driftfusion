%% Perform IMVS Measurement
par = pc('Input_files/3_layer_test_symmetric.csv');

soleq = equilibrate(par);

sol_IMVS = do_IMVS(soleq.no_ion, 1, 0.2, 4, 1, 200);

%% Plot the Voc as a function of t
dfplot.Voct(sol_IMVS);

%% Plot the generation rate as a function of x and t
x_ihalf = getvarihalf(sol_IMVS.x);

figure(200)
surf(x_ihalf, sol_IMVS.t, sol_IMVS.g)
xlabel('Position [cm]')
ylabel('Time [s]')
zlabel('Generation rate [cm-3s-1]')