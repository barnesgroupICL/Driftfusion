% Plot recombination
dfplot.Ux(JV.pin_srh_i.ill.r, 194)
hold on
dfplot.Ux(JV.hom_srh_i.ill.r, 170)
hold on
dfplot.Ux(JV.ntype_srh_i.ill.r, 170)
hold on
dfplot.Ux(JV.mobile_dope_i.ill.r, 170)
legend('intrinsic', 'homojunction', 'n-type', 'mobile dopants')
hold off
