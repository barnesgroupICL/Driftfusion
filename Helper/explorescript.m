
par.ptpd = pc('input_files/ptpd.csv');
ptpd_parex_dactive_light_ = explore.explore2par(par.ptpd, {'dcell(1,4)','Int'}, {[10e-7, 40e-7, 90e-7, 140e-7, 340e-7], logspace(-4,1,6)}, 400);