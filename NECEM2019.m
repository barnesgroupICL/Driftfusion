%% NECEM 2019 script
%  par.tio2 = pc('input_files/spiro_mapi_tio2_nosr.csv');
%  soleq.tio2 = equilibrate(par.tio2);
%  JV.tio2 = doJV(soleq.tio2.ion, 100e-3, 201, 1, 1, 0, 1.2, 2);
%  JV.tio2_no_ion = doJV(soleq.tio2.no_ion, 100e-3, 201, 1, 1, 0, 1.2, 2);
% 
%  dfplot.ELx_single(JV.tio2_no_ion.ill.f, 0);
%  dfplot.npx(JV.tio2_no_ion.ill.f, 0);

 par.pn = pc('input_files/pn_junction.csv');
 soleq.pn = equilibrate(par.pn);
%  JV.pn = doJV(soleq.pn.ion, 100e-3, 201, 1, 1, 0, 1.2, 2);
JV.pn_no_ion = doJV(soleq.pn.no_ion, 100e-3, 100, 1, 1, 0, 0.6, 1);

pn_dk_f = moviemake(JV.pn_no_ion.dk.f, @dfplot.ELnpx, 0, [1e4, 1e17], 'pn_dk_f_mov');

%  dfplot.ELx_single(JV.pn_no_ion.ill.f, 0);
%  dfplot.npx(JV.pn_no_ion.ill.f, 0);
