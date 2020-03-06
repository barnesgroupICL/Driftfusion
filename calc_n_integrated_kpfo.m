%% electrons
[kpfo_n.intrinsic, n_integrated.intrinsicl] = dfana.calckpfo_n(OC.intrinsic, 21e-7, 561e-7);
[kpfo_n.intrinsic_i, n_integrated.intrinsic_i] = dfana.calckpfo_n(OC.intrinsic_i, 21e-7, 561e-7);

[kpfo_n.hom, n_integrated.hom] = dfana.calckpfo_n(OC.hom, 21e-7, 561e-7);
[kpfo_n.hom_i, n_integrated.hom_i] = dfana.calckpfo_n(OC.hom_i, 21e-7, 561e-7);

[kpfo_n.ntype, n_integrated.ntype] = dfana.calckpfo_n(OC.ntype, 21e-7, 561e-7);
[kpfo_n.ntype_i, n_integrated.ntype_i] = dfana.calckpfo_n(OC.ntype_i, 21e-7, 561e-7);

[kpfo_n.mobile_dope, n_integrated.mobile_dope] = dfana.calckpfo_n(OC.mobile_dope, 21e-7, 561e-7);
[kpfo_n.mobile_dope_i, n_integrated.mobile_dope_i] = dfana.calckpfo_n(OC.mobile_dope_i, 21e-7, 561e-7);

%% holes
[kpfo_p.intrinsic, p_integrated.intrinsic] = dfana.calckpfo_p(OC.intrinsic, 21e-7, 561e-7);
[kpfo_p.intrinsic_i, p_integrated.intrinsic_i] = dfana.calckpfo_p(OC.intrinsic_i, 21e-7, 561e-7);

[kpfo_p.hom, p_integrated.hom] = dfana.calckpfo_p(OC.hom, 21e-7, 561e-7);
[kpfo_p.hom_i, p_integrated.hom_i] = dfana.calckpfo_p(OC.hom_i, 21e-7, 561e-7);

[kpfo_p.ntype, p_integrated.ntype] = dfana.calckpfo_p(OC.ntype, 21e-7, 561e-7);
[kpfo_p.ntype_i, p_integrated.ntype_i] = dfana.calckpfo_p(OC.ntype_i, 21e-7, 561e-7);

[kpfo_p.mobile_dope, p_integrated.mobile_dope] = dfana.calckpfo_p(OC.mobile_dope, 21e-7, 561e-7);
[kpfo_p.mobile_dope_i, p_integrated.mobile_dope_i] = dfana.calckpfo_p(OC.mobile_dope_i, 21e-7, 561e-7);

%% QFL splitting calculated from average active layer carrier densities
[deltaQFL.intrinsic, nbarpbar.intrinsic] = calcdeltaQFL(n_integrated.intrinsic, p_integrated.intrinsic, 2.50e20, 2.50e20, 540e-7, 1.5);
[deltaQFL.intrinsic_i, nbarpbar.intrinsic_i] = calcdeltaQFL(n_integrated.intrinsic_i, p_integrated.intrinsic_i, 2.50e20, 2.50e20, 540e-7, 1.5);

[deltaQFL.hom, nbarpbar.hom] = calcdeltaQFL(n_integrated.hom, p_integrated.hom, 2.50e20, 2.50e20, 540e-7, 1.5);
[deltaQFL.hom_i, nbarpbar.hom_i] = calcdeltaQFL(n_integrated.hom_i, p_integrated.hom_i, 2.50e20, 2.50e20, 540e-7, 1.5);

[deltaQFL.ntype, nbarpbar.ntype] = calcdeltaQFL(n_integrated.ntype, p_integrated.ntype, 2.50e20, 2.50e20, 540e-7, 1.5);
[deltaQFL.ntype_i, nbarpbar.ntype_i] = calcdeltaQFL(n_integrated.ntype_i, p_integrated.ntype_i, 2.50e20, 2.50e20, 540e-7, 1.5);

[deltaQFL.mobile_dope, nbarpbar.mobile_dope]  = calcdeltaQFL(n_integrated.mobile_dope, p_integrated.mobile_dope, 2.50e20, 2.50e20, 540e-7, 1.5);
[deltaQFL.mobile_dope_i, nbarpbar.mobile_dope_i] = calcdeltaQFL(n_integrated.mobile_dope_i, p_integrated.mobile_dope_i, 2.50e20, 2.50e20, 540e-7, 1.5);

%% Simulation Voc readout
VQFL.intrinsic = dfana.calcVQFL(OC.intrinsic);
VQFL.intrinsic = abs(VQFL.intrinsic(end));
VQFL.intrinsic_Nion1e19_i = dfana.calcVQFL(OC.intrinsic_i);
VQFL.intrinsic_Nion1e19_i = abs(VQFL.intrinsic_Nion1e19_i(end));
VQFL.intrinsic_Nion1e17_i = dfana.calcVQFL(OC.intrinsic_Nion1e17_i);
VQFL.intrinsic_Nion1e17_i = abs(VQFL.intrinsic_Nion1e17_i(end));

VQFL.hom = dfana.calcVQFL(OC.hom);
VQFL.hom = abs(VQFL.hom(end));
VQFL.hom_Nion1e19_i = dfana.calcVQFL(OC.hom_i);
VQFL.hom_Nion1e19_i = abs(VQFL.hom_Nion1e19_i(end));
VQFL.hom_Nion1e17_i = dfana.calcVQFL(OC.hom_Nion1e17_i);
VQFL.hom_Nion1e17_i = abs(VQFL.hom_Nion1e17_i(end));

VQFL.ntype = dfana.calcVQFL(OC.ntype);
VQFL.ntype = abs(VQFL.ntype(end));
VQFL.ntype_Nion1e19_i = dfana.calcVQFL(OC.ntype_i);
VQFL.ntype_Nion1e19_i = abs(VQFL.ntype_Nion1e19_i(end));
VQFL.ntype_Nion1e17_i = dfana.calcVQFL(OC.ntype_Nion1e17_i);
VQFL.ntype_Nion1e17_i = abs(VQFL.ntype_Nion1e17_i(end));

VQFL.mobile_dope = dfana.calcVQFL(OC.mobile_dope);
VQFL.mobile_dope = abs(VQFL.mobile_dope(end));
VQFL.mobile_dope_Nion1e19_i = dfana.calcVQFL(OC.mobile_dope_i);
VQFL.mobile_dope_Nion1e19_i = abs(VQFL.mobile_dope_Nion1e19_i(end));
