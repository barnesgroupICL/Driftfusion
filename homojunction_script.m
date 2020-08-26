%% Intrinsic, Nion= 1e19 cm-3
%% Intrinsic
par.intrinsic = pc('input_files/intrinsic.csv');
soleq.intrinsic = equilibrate(par.intrinsic);
JV.intrinsic = doJV(soleq.intrinsic.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.intrinsic_i = doJV(soleq.intrinsic.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
OC.intrinsic = lightonRs(soleq.intrinsic.el, 1, 1, 0, 1e6, 200);
OC.intrinsic_i = lightonRs(soleq.intrinsic.ion, 1, -1, 1, 1e6, 200);

%% n-type
par.ntype = pc('input_files/n_type.csv');
soleq.ntype = equilibrate(par.ntype);
JV.ntype = doJV(soleq.ntype.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.ntype_i = doJV(soleq.ntype.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
OC.ntype = lightonRs(soleq.ntype.el, 1, 1, 0, 1e6, 200);
OC.ntype_i = lightonRs(soleq.ntype.ion, 1, -1, 0, 1e6, 200);

%% Homojunction
par.hom = pc('input_files/homojunction.csv');
soleq.hom = equilibrate(par.hom);
JV.hom = doJV(soleq.hom.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.hom_i = doJV(soleq.hom.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
OC.hom = lightonRs(soleq.hom.el, 1, 1, 0, 1e6, 200);
OC.hom_i = lightonRs(soleq.hom.ion, 1, -1, 0, 1e6, 200);

%% Mobile dopants
par.mobile_dope = pc('input_files/mobile_dopants.csv');
soleq.mobile_dope = equilibrate(par.mobile_dope);
JV.mobile_dope = doJV(soleq.mobile_dope.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.mobile_dope_i = doJV(soleq.mobile_dope.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
OC.mobile_dope = lightonRs(soleq.mobile_dope.el, 1, 1, 0, 1e6, 200);
OC.mobile_dope_i = lightonRs(soleq.mobile_dope.ion, 1, -1, 0, 1e6, 200);

%% Intrinsic, Nion= 1e17 cm-3
par.intrinsic_Nion1e17 = pc('input_files/intrinsic.csv');
par.intrinsic_Nion1e17.Ncat = [1e17,1e17,1e17,1e17,1e17];
par.intrinsic_Nion1e17.Nani = [1e17,1e17,1e17,1e17,1e17];
par.intrinsic_Nion1e17 = refresh_device(par.intrinsic_Nion1e17);

soleq.intrinsic_Nion1e17 = equilibrate(par.intrinsic_Nion1e17);
% JV.intrinsic_Nion1e17 = doJV(soleq.intrinsic_Nion1e17.el, 1e-3, 100, 1, 0, 0, 1.2, 2);
JV.intrinsic_Nion1e17_i = doJV(soleq.intrinsic_Nion1e17.ion, 1e-3, 100, 1, 1, 0, 1.2, 2);
% OC.intrinsic_Nion1e17 = lightonRs(soleq.intrinsic_Nion1e17.el, 1, 1, 0, 1e6, 200);
OC.intrinsic_Nion1e17_i = lightonRs(soleq.intrinsic_Nion1e17.ion, 1, -1, 1, 1e6, 200);

%% Homojunction, Nion= 1e17 cm-3
par.hom_Nion1e17 = pc('input_files/homojunction.csv');
par.hom_Nion1e17.Ncat = [1e17,1e17,1e17,1e17,1e17,1e17,1e17];
par.hom_Nion1e17.Nani = [1e17,1e17,1e17,1e17,1e17,1e17,1e17];
par.hom_Nion1e17 = refresh_device(par.hom_Nion1e17);

soleq.hom_Nion1e17 = equilibrate(par.hom_Nion1e17);
% JV.hom_Nion1e17 = doJV(soleq.hom_Nion1e17.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.hom_Nion1e17_i = doJV(soleq.hom_Nion1e17.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
% OC.hom_Nion1e17 = lightonRs(soleq.hom_Nion1e17.el, 1, 1, 0, 1e6, 200);
OC.hom_Nion1e17_i = lightonRs(soleq.hom_Nion1e17.ion, 1, -1, 1, 1e6, 200);

%% n-type, Nion = 1e17 cm-3
par.ntype_Nion1e17 = pc('input_files/n_type.csv');
par.ntype_Nion1e17.Ncat = [1e17,1e17,1e17,1e17,1e17];
par.ntype_Nion1e17.Nani = [1e17,1e17,1e17,1e17,1e17];
par.ntype_Nion1e17 = refresh_device(par.ntype_Nion1e17);

soleq.ntype_Nion1e17 = equilibrate(par.ntype_Nion1e17);
% JV.ntype_Nion1e17 = doJV(soleq.ntype_Nion1e17.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.ntype_Nion1e17_i = doJV(soleq.ntype_Nion1e17.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
% OC.ntype_Nion1e17 = lightonRs(soleq.ntype_Nion1e17.el, 1, 1, 0, 1e6, 200);
OC.ntype_Nion1e17_i = lightonRs(soleq.ntype_Nion1e17.ion, 1, -1, 1, 1e6, 200);

%% Intrinsic, Nion= 1e18 cm-3
par.intrinsic_Nion1e18 = pc('input_files/intrinsic.csv');
par.intrinsic_Nion1e18.Ncat = [1e18,1e18,1e18,1e18,1e18];
par.intrinsic_Nion1e18.Nani = [1e18,1e18,1e18,1e18,1e18];
par.intrinsic_Nion1e18 = refresh_device(par.intrinsic_Nion1e18);

soleq.intrinsic_Nion1e18 = equilibrate(par.intrinsic_Nion1e18);
% JV.intrinsic_Nion1e18 = doJV(soleq.intrinsic_Nion1e18.el, 1e-3, 100, 1, 0, 0, 1.2, 2);
JV.intrinsic_Nion1e18_i = doJV(soleq.intrinsic_Nion1e18.ion, 1e-3, 100, 1, 1, 0, 1.2, 2);
% OC.intrinsic_Nion1e18 = lightonRs(soleq.intrinsic_Nion1e18.el, 1, 1, 0, 1e6, 200);
OC.intrinsic_Nion1e18_i = lightonRs(soleq.intrinsic_Nion1e18.ion, 1, -1, 1, 1e6, 200);

%% Homojunction, Nion= 1e18 cm-3
par.hom_Nion1e18 = pc('input_files/homojunction.csv');
par.hom_Nion1e18.Ncat = [1e18,1e18,1e18,1e18,1e18,1e18,1e18];
par.hom_Nion1e18.Nani = [1e18,1e18,1e18,1e18,1e18,1e18,1e18];
par.hom_Nion1e18 = refresh_device(par.hom_Nion1e18);

soleq.hom_Nion1e18 = equilibrate(par.hom_Nion1e18);
% JV.hom_Nion1e18 = doJV(soleq.hom_Nion1e18.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.hom_Nion1e18_i = doJV(soleq.hom_Nion1e18.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
% OC.hom_Nion1e18 = lightonRs(soleq.hom_Nion1e18.el, 1, 1, 0, 1e6, 200);
OC.hom_Nion1e18_i = lightonRs(soleq.hom_Nion1e18.ion, 1, -1, 1, 1e6, 200);

%% n-type, Nion = 1e18 cm-3
par.ntype_Nion1e18 = pc('input_files/n_type.csv');
par.ntype_Nion1e18.Ncat = [1e18,1e18,1e18,1e18,1e18];
par.ntype_Nion1e18.Nani = [1e18,1e18,1e18,1e18,1e18];
par.ntype_Nion1e18 = refresh_device(par.ntype_Nion1e18);

soleq.ntype_Nion1e18 = equilibrate(par.ntype_Nion1e18);
% JV.ntype_Nion1e18 = doJV(soleq.ntype_Nion1e18.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.ntype_Nion1e18_i = doJV(soleq.ntype_Nion1e18.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);
% OC.ntype_Nion1e18 = lightonRs(soleq.ntype_Nion1e18.el, 1, 1, 0, 1e6, 200);
OC.ntype_Nion1e18_i = lightonRs(soleq.ntype_Nion1e18.ion, 1, -1, 1, 1e6, 200);

%% J turn ons, Nion = 1e19 cm-3
sol_relax.intrinsic_el = jumptoV(soleq.intrinsic.el, -0.85, 10, 0, 1, 0, 0);
sol_relax.intrinsic_ion = jumptoV(soleq.intrinsic.ion, -0.92, 10, 1, 1, 0, 0);
sol_relax.hom_el = jumptoV(soleq.hom.el, -0.93, 10, 0, 1, 0, 0);
sol_relax.hom_ion = jumptoV(soleq.hom.ion, -0.95, 10, 1, 1, 0, 0);

%% Steady state increasing voltage
sol_relax.intrinsic_Nion1e18_el_m1V = jumptoV(soleq.intrinsic_Nion1e18.el, 1.0, 0.1, 0, 0, 0, 0);
sol_relax.intrinsic_Nion1e18_el_m0p5V = jumptoV(soleq.intrinsic_Nion1e18.el, 0.5, 0.1, 0, 0, 0, 0);
sol_relax.intrinsic_Nion1e18_el_0p25V = jumptoV(soleq.intrinsic_Nion1e18.el, -0.25, 0.1, 0, 0, 0, 0);
sol_relax.intrinsic_Nion1e18_el_0p5V = jumptoV(soleq.intrinsic_Nion1e18.el, -0.5, 0.1, 0, 0, 0, 0);

sol_relax.intrinsic_Nion1e18_ion_m1V = jumptoV(soleq.intrinsic_Nion1e18.ion, 1.0, 10, 1, 0, 0, 0);
sol_relax.intrinsic_Nion1e18_ion_m0p5V = jumptoV(soleq.intrinsic_Nion1e18.ion, 0.5, 10, 1, 0, 0, 0);
sol_relax.intrinsic_Nion1e18_ion_0p25V = jumptoV(soleq.intrinsic_Nion1e18.ion, -0.25, 10, 1, 0, 0, 0);
sol_relax.intrinsic_Nion1e18_ion_0p5V = jumptoV(soleq.intrinsic_Nion1e18.ion, -0.5, 10, 1, 0, 0, 0);

sol_relax.hom_Nion1e18_el_m1V = jumptoV(soleq.hom_Nion1e18.el, 1.0, 0.1, 0, 0, 0, 0);
sol_relax.hom_Nion1e18_el_m0p5V  = jumptoV(soleq.hom_Nion1e18.el, 0.5, 0.1, 0, 0, 0, 0);
sol_relax.hom_Nion1e18_el_0p25V = jumptoV(soleq.hom_Nion1e18.el, -0.25, 0.1, 0, 0, 0, 0);
sol_relax.hom_Nion1e18_el_0p5V = jumptoV(soleq.hom_Nion1e18.el, -0.5, 0.1, 0, 0, 0, 0);

sol_relax.hom_Nion1e18_ion_m1V = jumptoV(soleq.hom_Nion1e18.ion, 1.0, 10, 1, 0, 0, 0);
sol_relax.hom_Nion1e18_ion_m0p5V = jumptoV(soleq.hom_Nion1e18.ion, 0.5, 10, 1, 0, 0, 0);
sol_relax.hom_Nion1e18_ion_0p25V = jumptoV(soleq.hom_Nion1e18.ion, -0.25, 10, 1, 0, 0, 0);
sol_relax.hom_Nion1e18_ion_0p5V = jumptoV(soleq.hom_Nion1e18.ion, -0.5, 10, 1, 0, 0, 0);


