%% Intrinsic, Nion = 1e18 cm-3
dfplot.Fx(sol_relax.intrinsic_Nion1e18_el_m1V)
hold on
dfplot.Fx(sol_relax.intrinsic_Nion1e18_el_m0p5V)
hold on
dfplot.Fx(soleq.intrinsic_Nion1e18.ion)
hold on
dfplot.Fx(sol_relax.intrinsic_Nion1e18_el_0p25V)
hold on
dfplot.Fx(sol_relax.intrinsic_Nion1e18_el_0p5V)
legend('HTL','int1','Active','int2','ETL','Intrinsic, Vapp = -1 V','Intrinsic, Vapp = -0.5 V',...
    'Intrinsic, Vapp = 0 V','Intrinsic, Vapp = 0.25 V','Intrinsic, Vapp = 0.5 V')
hold off
ylim([-0.5e5,2e5])

%% Homojunction, Nion = 1e18 cm-3
dfplot.Fx(sol_relax.hom_Nion1e18_el_m1V)
hold on
dfplot.Fx(sol_relax.hom_Nion1e18_el_m0p5V)
hold on
dfplot.Fx(soleq.hom_Nion1e18.el)
hold on
dfplot.Fx(sol_relax.hom_Nion1e18_el_0p25V)
hold on
dfplot.Fx(sol_relax.hom_Nion1e18_el_0p5V)
legend('HTL','int1','n-type','int2','p-type','int3','ETL','Homojunction, Vapp = -1 V','Homojunction, Vapp = -0.5 V',...
    'Homojunction, Vapp = 0 V','Homojunction, Vapp = 0.25 V','Homojunction, Vapp = 0.5 V')
hold off
ylim([-0.5e5,2e5])

%% Homojunction, Nion = 1e18 cm-3
dfplot.Fx(sol_relax.hom_Nion1e18_ion_m1V)
hold on
dfplot.Fx(sol_relax.hom_Nion1e18_ion_m0p5V)
hold on
dfplot.Fx(soleq.hom_Nion1e18.ion)
hold on
dfplot.Fx(sol_relax.hom_Nion1e18_ion_0p25V)
hold on
dfplot.Fx(sol_relax.hom_Nion1e18_ion_0p5V)
legend('HTL','int1','n-type','int2','p-type','int3','ETL','Homojunction, Vapp = -1 V','Homojunction, Vapp = -0.5 V',...
    'Homojunction, Vapp = 0 V','Homojunction, Vapp = 0.25 V','Homojunction, Vapp = 0.5 V')
hold off
ylim([-0.5e5,2e5])
