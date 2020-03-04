%% QFL splitting calculated from average active layer carrier densities
[deltaQFL_midpoint.intrinsic, np_midpoint.intrinsic] = calcdeltaQFL_midpoint(OC.intrinsic, 2.50e20, 2.50e20, 540e-7, 1.5);
[deltaQFL_midpoint.intrinsic_i, np_midpoint.intrinsic_i] = calcdeltaQFL_midpoint(OC.intrinsic_i, 2.50e20, 2.50e20, 540e-7, 1.5);

[deltaQFL_midpoint.hom, np_midpoint.hom] = calcdeltaQFL_midpoint(OC.hom,  2.50e20, 2.50e20, 540e-7, 1.5);
[deltaQFL_midpoint.hom_i, np_midpoint.hom_i] = calcdeltaQFL_midpoint(OC.hom_i, 2.50e20, 2.50e20, 540e-7, 1.5);

[deltaQFL_midpoint.ntype, np_midpoint.ntype] = calcdeltaQFL_midpoint(OC.ntype,  2.50e20, 2.50e20, 540e-7, 1.5);
[deltaQFL_midpoint.ntype_i, np_midpoint.ntype_i] = calcdeltaQFL_midpoint(OC.ntype_i,  2.50e20, 2.50e20, 540e-7, 1.5);

[deltaQFL_midpoint.mobile_dope, np_midpoint.mobile_dope]  = calcdeltaQFL_midpoint(OC.mobile_dope, 2.50e20, 2.50e20, 540e-7, 1.5);
[deltaQFL_midpoint.mobile_dope_i, np_midpoint.mobile_dope_i] = calcdeltaQFL_midpoint(OC.mobile_dope_i, 2.50e20, 2.50e20, 540e-7, 1.5);

