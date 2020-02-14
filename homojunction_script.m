%% Intrinsic
par.intrinsic = pc('input_files/intrinsic.csv');
soleq.intrinsic = equilibrate(par.intrinsic);
JV.intrinsic = doJV(soleq.intrinsic.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.intrinsic_i = doJV(soleq.intrinsic.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);

%% n-type
par.ntype = pc('input_files/n_type.csv');
soleq.ntype = equilibrate(par.ntype);
JV.ntype = doJV(soleq.ntype.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.ntype_i = doJV(soleq.ntype.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);

%% Homojunction
par.hom = pc('input_files/homojunction.csv');
soleq.hom = equilibrate(par.hom);
JV.hom = doJV(soleq.hom.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.hom_i = doJV(soleq.hom.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);

%% Mobile dopants
par.mobile_dope = pc('input_files/mobile_dopants.csv');
soleq.mobile_dope = equilibrate(par.mobile_dope);
JV.mobile_dope = doJV(soleq.mobile_dope.el, 1e-3, 100, 1, 0, 0, -1.2, 2);
JV.mobile_dope_i = doJV(soleq.mobile_dope.ion, 1e-3, 100, 1, 1, 0, -1.2, 2);

%%
par.intrinsic_Nion1e17 = pc('input_files/Nion1e17 intrinsic.csv');
soleq.intrinsic_Nion1e17 = equilibrate(par.intrinsic_Nion1e17);
JV.intrinsic_Nion1e17 = doJV(soleq.intrinsic_Nion1e17.el, 1e-3, 100, 1, 0, 0, 1.2, 2);
JV.intrinsic_Nion1e17_i = doJV(soleq.intrinsic_Nion1e17.ion, 1e-3, 100, 1, 1, 0, 1.2, 2);

%%
par.hom_Nion1e17 = pc('input_files/Nion1e17 homojunction.csv');
soleq.hom_Nion1e17 = equilibrate(par.hom_Nion1e17);
JV.hom_Nion1e17 = doJV(soleq.hom_Nion1e17.el, 1e-3, 100, 1, 0, 0, 1.2, 2);
JV.hom_Nion1e17_i = doJV(soleq.hom_Nion1e17.ion, 1e-3, 100, 1, 1, 0, 1.2, 2);

%%
par.ntype_Nion1e17 = pc('input_files/Nion1e17 ntype.csv');
soleq.ntype_Nion1e17 = equilibrate(par.ntype_Nion1e17);
JV.ntype_Nion1e17 = doJV(soleq.ntype_Nion1e17.el, 1e-3, 100, 1, 0, 0, 1.2, 2);
JV.ntype_Nion1e17_i = doJV(soleq.ntype_Nion1e17.ion, 1e-3, 100, 1, 1, 0, 1.2, 2);

%%
