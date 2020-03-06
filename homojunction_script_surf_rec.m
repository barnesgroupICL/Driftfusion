%% Intrinsic
par.intrinsic = pc('input_files/intrinsic_surf_rec.csv');
soleq.intrinsic = equilibrate(par.intrinsic);
JV.intrinsic = doJV(soleq.intrinsic.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.intrinsic_i = doJV(soleq.intrinsic.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
OC.intrinsic = lightonRs(soleq.intrinsic.el, 1, 1, 0, 1e6, 200);
OC.intrinsic_i = lightonRs(soleq.intrinsic.ion, 1, -1, 1, 1e6, 200);

%% Intrinsic, Nion= 1e17 cm-3
par.intrinsic_Nion1e17 = pc('input_files/intrinsic_surf_rec.csv');
par.intrinsic_Nion1e17.Ncat = [1e17,1e17,1e17,1e17,1e17];
par.intrinsic_Nion1e17.Nani = [1e17,1e17,1e17,1e17,1e17];
par.intrinsic_Nion1e17 = refresh_device(par.intrinsic_Nion1e17);

soleq.intrinsic_Nion1e17 = equilibrate(par.intrinsic_Nion1e17);
JV.intrinsic_Nion1e17 = doJV(soleq.intrinsic_Nion1e17.el, 1e-3, 100, 1, 0, 0, 1.2, 2);
JV.intrinsic_Nion1e17_i = doJV(soleq.intrinsic_Nion1e17.ion, 1e-3, 100, 1, 1, 0, 1.2, 2);
OC.intrinsic_Nion1e17 = lightonRs(soleq.intrinsic_Nion1e17.el, 1, 1, 0, 1e6, 200);
OC.intrinsic_Nion1e17_i = lightonRs(soleq.intrinsic_Nion1e17.ion, 1, -1, 1, 1e6, 200);

