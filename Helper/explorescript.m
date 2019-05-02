par.ptpd = pc('input_files/ptpd.csv');
ptpd_parex_dactive_light_ = explore.explore2par(par.ptpd, {'dcell(1,4)','Int'}, {[10e-7, 40e-7, 90e-7, 140e-7, 340e-7, 740e-7], logspace(-3,1,6)}, 400);

par.pcbm = pc('input_files/pcbm.csv');
pcbm_parex_dactive_light_ = explore.explore2par(par.pcbm, {'dcell(1,4)','Int'}, {[10e-7, 40e-7, 90e-7, 140e-7, 340e-7, 740e-7], logspace(-3,1,6)}, 400);

par.ptpd_Na19 = pc('input_files/ptpd.csv');
par.ptpd_Na19.Nion = 1e19;

ptpd_parex_dactive_light_ = explore.explore2par(par.ptpd_Na19, {'dcell(1,4)','Int'}, {[10e-7, 40e-7, 90e-7, 140e-7, 340e-7, 740e-7], logspace(-3,1,6)}, 400);

par.pcbm_Na19 = pc('input_files/pcbm.csv');
par.pcbm_Na19.Nion = 1e19;
pcbm_parex_dactive_light_ = explore.explore2par(par.pcbm_Na19, {'dcell(1,4)','Int'}, {[10e-7, 40e-7, 90e-7, 140e-7, 340e-7, 740e-7], logspace(-3,1,6)}, 400);