
par.p = pc('input_files/tio2.csv');
tio2_parex_dactive_light_ = explore.explore2par(par.tio2, {'dcell(1,5)','Int'}, {[10e-7, 40e-7, 90e-7, 140e-7, 340e-7], logspace(-4,1,6)}, 400);