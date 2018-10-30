% parexsol_NA_taun_PhiA_ptpd_NI0 = explore.explore2par(par.ptpd.base, {'taun(1)', 'E0(1)'}, {logspace(-13, -6, 8); -5.0:0.1:-4.3});
% parexsol_NA_taun_PhiA_pedot_NI0 = explore.explore2par(par.pedot.base, {'taun(1)', 'E0(1)'}, {logspace(-13, -6, 8); -5.0:0.1:-4.3})
% parexsol_NA_NI_ptpd = explore.explore2par(par.ptpd.base, {'NI', 'E0(1)'}, {logspace(16, 20, 5); -5.0:0.1:-4.7});
% parexsol_NA_NI_pedot = explore.explore2par(par.pedot.base, {'NI', 'E0(1)'}, {logspace(16, 20, 5); -5.0:0.1:-4.7});

parexsol_PhiA_NI_ptpd = explore.explore2par(par.ptpd.base, {'taun(1)', 'PhiA'}, {logspace(-13, -6, 8); -5.0:0.1:-4.3});
parexsol_PhiA_NI_pedot = explore.explore2par(par.pedot.base, {'taun(1)', 'PhiA'}, {logspace(-13, -6, 8); -5.0:0.1:-4.3});