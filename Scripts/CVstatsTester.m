par = pc('Input_files/1_layer_test.csv');

sol = equilibrate(par);

sol_CV1 = doCV(sol.el, 1, 0.7, 0.9, 0, 100e-3, 3, 241);
sol_CV2 = doCV(sol.el, 1, 0, 0.9, 0, 100e-3, 3, 241);
sol_CV3 = doCV(sol.el, 1, 0, 0.9, 0, 100e-3, 1, 241);