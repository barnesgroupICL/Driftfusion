par = pc('Input_files/1_layer_test.csv');

sol = equilibrate(par);

sol_CV1 = doCV(sol.el, 1, 0.7, 1.0, 0, 100e-3, 3, 241);
sol_CV2 = doCV(sol.el, 1, 0, 1.0, 0, 100e-3, 3, 241);
sol_CV3 = doCV(sol.el, 1, 0, 1.0, 0, 100e-3, 1, 241);
sol_CV4 = doCV(sol.el, 1, 0, 0.9, 0, 100e-3, 3, 241);

stats1 = CVstats(sol_CV1)
stats2 = CVstats(sol_CV2)
stats3 = CVstats(sol_CV3)
stats4 = CVstats(sol_CV4)