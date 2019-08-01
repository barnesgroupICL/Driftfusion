%% NECEM 2019 script
%  par.tio2 = pc('input_files/spiro_mapi_tio2_nosr.csv');
%  soleq.tio2 = equilibrate(par.tio2);
%  JV.tio2 = doJV(soleq.tio2.ion, 100e-3, 201, 1, 1, 0, 1.2, 2);
%  JV.tio2_no_ion = doJV(soleq.tio2.no_ion, 100e-3, 201, 1, 1, 0, 1.2, 2);
% JV.tio2_ionic_branch = doJV(soleq.tio2.ion, 100e-3, 201, 1, 1, 0, 0.8, 1);
% JV.tio2_slowscan = doJV(soleq.tio2.ion, 1e-3, 101, 1, 1e10, 0, 1.2, 2);
% %  dfplot.ELx_single(JV.tio2_no_ion.ill.f, 0);
% %  dfplot.npx(JV.tio2_no_ion.ill.f, 0);
% 
 par.pn = pc('input_files/pn_junction.csv');
 soleq.pn = equilibrate(par.pn);
%  JV.pn = doJV(soleq.pn.ion, 100e-3, 201, 1, 1, 0, 1.2, 2);
JV.pn_no_ion = doJV(soleq.pn.no_ion, 100e-3, 100, 1, 1, 0, 0.6, 1);
% 
% pn_dk_f = moviemake(JV.pn_no_ion.dk.f, @dfplot.ELnpx, 0, [1e4, 1e17], 'pn_dk_f_mov');

%  dfplot.ELx_single(JV.pn_no_ion.ill.f, 0);
%  dfplot.npx(JV.pn_no_ion.ill.f, 0);

% par.pin = pc('input_files/Evidence_for_ion_migration_SRH.csv');
% soleq.pin = equilibrate(par.pin);
%  JV.pn = doJV(soleq.pn.ion, 100e-3, 201, 1, 1, 0, 1.2, 2);
% JV.pin_no_ion_0p1Vs = doJV(soleq.pin.ion, 0.1, 100, 1, 1, 0, 1, 2);
% JV.pin_no_ion_1Vs = doJV(soleq.pin.ion, 1, 100, 1, 1, 0, 1, 2);
% JV.pin_no_ion_10Vs = doJV(soleq.pin.ion, 10, 100, 1, 1, 0, 1, 2);
%JV.pin_no_ion_ss = doJV(soleq.pin.ion, 1, 100, 1, 1e10, 0, 1, 2);

% par.pn_ion = pc('input_files/pn_junction_ion.csv');
% soleq.pn_ion = equilibrate(par.pn_ion);
% %  JV.pn_ion = doJV(soleq.pn_ion.ion, 100e-3, 201, 1, 1, 0, 1.2, 2);
% JV.pn_ion = doJV(soleq.pn_ion.ion, 100e-3, 100, 1, 1, 0, 0.6, 1);
% pn_ion_dk_f = moviemake(JV.pn_ion.dk.f, @dfplot.ELnpx, 0, [1e4, 1e17], 'pn_ion_dk_f_mov');