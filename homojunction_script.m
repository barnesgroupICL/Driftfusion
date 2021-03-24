%% Intrinsic, Nion= 1e19 cm-3
%% Intrinsic
par_intrinsic = pc('input_files/intrinsic.csv');
soleq_intrinsic = equilibrate(par_intrinsic);
JV_intrinsic = doJV(soleq_intrinsic.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV_intrinsic_i = doJV(soleq_intrinsic.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
OC_intrinsic = lightonRs(soleq_intrinsic.el, 1, 1, 0, 1e6, 200);
OC_intrinsic_i = lightonRs(soleq_intrinsic.ion, 1, -100, 1, 1e6, 200);

%% n-type
par_ntype = pc('input_files/n_type.csv');
soleq_ntype = equilibrate(par_ntype);
JV_ntype = doJV(soleq_ntype.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV_ntype_i = doJV(soleq_ntype.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
OC_ntype = lightonRs(soleq_ntype.el, 1, 1, 0, 1e6, 200);
OC_ntype_i = lightonRs(soleq_ntype.ion, 1, -100, 0, 1e6, 200);

%% Homojunction
par_hom = pc('input_files/homojunction.csv');
soleq_hom = equilibrate(par_hom);
JV_hom = doJV(soleq_hom.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV_hom_i = doJV(soleq_hom.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
OC_hom = lightonRs(soleq_hom.el, 1, 1, 0, 1e6, 200);
OC_hom_i = lightonRs(soleq_hom.ion, 1, -100, 0, 1e6, 200);

%% Homojunction, Nion= 1e17 cm-3
par_hom_Nion1e17 = pc('input_files/homojunction.csv');
par_hom_Nion1e17.Ncat = [1e17,1e17,1e17,1e17,1e17,1e17,1e17];
par_hom_Nion1e17.Nani = [1e17,1e17,1e17,1e17,1e17,1e17,1e17];
par_hom_Nion1e17 = refresh_device(par_hom_Nion1e17);

soleq_hom_Nion1e17 = equilibrate(par_hom_Nion1e17);
JV_hom_Nion1e17_i = doJV(soleq_hom_Nion1e17.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
OC_hom_Nion1e17_i = lightonRs(soleq_hom_Nion1e17.ion, 1, -100, 1, 1e6, 200);

%% Homojunction, Nion= 1e18 cm-3
par_hom_Nion1e18 = pc('input_files/homojunction.csv');
par_hom_Nion1e18.Ncat = [1e18,1e18,1e18,1e18,1e18,1e18,1e18];
par_hom_Nion1e18.Nani = [1e18,1e18,1e18,1e18,1e18,1e18,1e18];
par_hom_Nion1e18 = refresh_device(par_hom_Nion1e18);

soleq_hom_Nion1e18 = equilibrate(par_hom_Nion1e18);
JV_hom_Nion1e18_i = doJV(soleq_hom_Nion1e18.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
OC_hom_Nion1e18_i = lightonRs(soleq_hom_Nion1e18.ion, 1, -100, 1, 1e6, 200);

%% J turn ons, Nion = 1e19 cm-3
sol_relax.intrinsic_el = jumptoV(soleq_intrinsic.el, -0.85, 10, 0, 1, 0, 0);
sol_relax.intrinsic_ion = jumptoV(soleq_intrinsic.ion, -0.92, 10, 1, 1, 0, 0);
sol_relax.hom_el = jumptoV(soleq_hom.el, -0.93, 10, 0, 1, 0, 0);
sol_relax.hom_ion = jumptoV(soleq_hom.ion, -0.95, 10, 1, 1, 0, 0);




