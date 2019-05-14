% par.ptpd_Na18 = pc('input_files/ptpd.csv');
% ptpd_Na18_parex_dactive_light_85 = explore.explore2par(par.ptpd_Na18, {'dcell(1,4)','Int'}, {[40e-7, 140e-7, 240e-7, 340e-7, 440e-7, 540e-7, 640e-7, 740e-7], logspace(-3,1,5)}, 400);

%par.pcbm = pc('input_files/pcbm.csv');
% pcbm_parex_dactive_light_65 = explore.explore2par(par.pcbm, {'dcell(1,4)','Int'}, {[10e-7, 40e-7, 90e-7, 140e-7, 340e-7, 740e-7], logspace(-3,1,5)}, 400);
 
% par.ptpd_Na19 = pc('input_files/ptpd.csv');
% par.ptpd_Na19.Nion = par.ptpd_Na19.Nion*10;
% par.ptpd_Na19.dev = pc.builddev(par.ptpd_Na19);
% ptpd_Na19_parex_dactive_light_85 = explore.explore2par(par.ptpd_Na19, {'dcell(1,4)','Int'}, {[40e-7, 140e-7, 240e-7, 340e-7, 440e-7, 540e-7, 640e-7, 740e-7], logspace(-3,1,5)}, 400);

% par.pcbm_Na19 = pc('input_files/pcbm.csv');
% par.pcbm_Na19.Nion = par.pcbm_Na19.Nion*10;
% par.pcbm_Na19.dev = pc.builddev(par.pcbm_Na19);
% pcbm_Na19_parex_dactive_light_65 = explore.explore2par(par.pcbm_Na19, {'dcell(1,4)','Int'}, {[10e-7, 40e-7, 90e-7, 140e-7, 340e-7, 740e-7], logspace(-3,1,5)}, 400);
% 
% par.ptpd_Na17 = pc('input_files/ptpd.csv');
% par.ptpd_Na17.Nion = par.ptpd_Na17.Nion/10;
% par.ptpd_Na17.dev = pc.builddev(par.ptpd_Na17);
% ptpd_Na17_parex_dactive_light_85 = explore.explore2par(par.ptpd_Na17, {'dcell(1,4)','Int'}, {[40e-7, 140e-7, 240e-7, 340e-7, 440e-7, 540e-7, 640e-7, 740e-7], logspace(-3,1,5)}, 400);

%% Equilibrium sols for charge extraction
par.ptpd_Na17 = pc('input_files/ptpd.csv');
par.ptpd_Na17.Nion = par.ptpd_Na17.Nion/10;
par.ptpd_Na17.dev = pc.builddev(par.ptpd_Na17);
eq_ptpd_Na17_dactive_light_85 = explore.explore2par(par.ptpd_Na17, {'dcell(1,4)','Int'}, {[40e-7, 140e-7, 240e-7, 340e-7, 440e-7, 540e-7, 640e-7, 740e-7], 0}, 400);

par.ptpd_Na18 = pc('input_files/ptpd.csv');
eq_ptpd_Na18_dactive_light_85 = explore.explore2par(par.ptpd_Na18, {'dcell(1,4)','Int'}, {[40e-7, 140e-7, 240e-7, 340e-7, 440e-7, 540e-7, 640e-7, 740e-7], 0}, 400);

par.ptpd_Na19 = pc('input_files/ptpd.csv');
par.ptpd_Na19.Nion = par.ptpd_Na19.Nion*10;
par.ptpd_Na19.dev = pc.builddev(par.ptpd_Na19);
eq_ptpd_Na19_dactive_light_85 = explore.explore2par(par.ptpd_Na19, {'dcell(1,4)','Int'}, {[40e-7, 140e-7, 240e-7, 340e-7, 440e-7, 540e-7, 640e-7, 740e-7], 0}, 400);


save('ws_explore.mat')  % Save workspace