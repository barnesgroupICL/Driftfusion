%% Perform IMVS Measurement
par = pc('Input_files/3_layer_test_symmetric.csv');

soleq = equilibrate(par);

sol_IMVS = do_IMVS(soleq.no_ion, 1, 0.2, 4, 1, 1000);

dfplot.Voct(sol_IMVS);

surf(sol_IMVS.x, sol_IMVS.t, sol_IMVS.g)
xlabel('Position [cm]')
ylabel('Time [s]')
zlabel('Generation rate [cm-3s-1]')