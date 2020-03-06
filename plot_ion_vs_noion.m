%% JVs
dfplot.JV(JV.intrinsic_i, 2)
figure(4)
ylim([-10e-3, 30e-3])
hold on
dfplot.JV(JV.ntype_i, 2)
figure(4)
hold on
dfplot.JV(JV.hom_i, 2)
figure(4)
hold on
dfplot.JV(JV.mobile_dope_i, 2)
figure(4)
hold on

%% Equilibrium, Mobile dopants, V np x
dfplot.Vnpx(soleq.mobile_dope.el)
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.Vnpx(soleq.mobile_dope.ion)
subplot(2, 1, 1)
hold off
legend('HTL','int1','p-type perovskite','int2','n-type perovskite','int3','ETL',...
        'static dopants', 'mobile dopants')
hold off
subplot(2, 1, 2)
hold off

%% Equilibrium, Mobile dopants, rhox Fx Vx
dfplot.rhoxFxVx(soleq.mobile_dope.el)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(soleq.mobile_dope.ion)
subplot(3, 1, 3)
hold off
legend('HTL','int1','p-type perovskite','int2','n-type perovskite','int3','ETL',...
        'static dopants', 'mobile dopants')
hold off
subplot(3, 1, 2)
hold off
subplot(3, 1, 3)
hold off

%% OC, Mobile dopants, rhox Fx Vx
dfplot.rhoxFxVx(OC.mobile_dope)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC.mobile_dope_i)
subplot(3, 1, 3)
hold off
legend('HTL','int1','Active','int2','ETL','static dopants','mobile dopants')
hold off
subplot(3, 1, 2)
hold off
subplot(3, 1, 3)
hold off

%% OC,EL np x
dfplot.ELnpx(OC.intrinsic)
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.ELnpx(OC.intrinsic_i)
subplot(2, 1, 1)
hold off
legend('HTL','int1','Active','int2','ETL','intrinsic \it{N}_{SD}=0cm^{-3}','intrinsic \it{N}_{SD}=1e19cm^{-3}')
subplot(2, 1, 2)
hold off

%% OC,ac np F x, intrinsic Nion = 0, 1e17, 1e19
dfplot.rhoacnpx(OC.intrinsic_i)
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.rhoacnpx(OC.intrinsic_Nion1e17_i)
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.rhoacnpx(OC.intrinsic)
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.rhoacnpx(OC.ntype)
subplot(2, 1, 1)
legend('HTL','int1','Active','int2','ETL','\it{N}_{SD} = 10^{19} cm^{-3}',...
    '\it{N}_{SD} = 10^{17} cm^{-3}, \it{p}','\it{N}_{SD} = 0 cm^{-3}, \it{p}', 'n-type')
ylim([-2e17, 4e17])
hold off
subplot(2, 1, 2)
legend('HTL','int1','Active','int2','ETL','\it{N}_{SD} = 10^{19} cm^{-3}, \it{n}','\it{N}_{SD} = 10^{19} cm^{-3}, \it{p}',...
    '\it{N}_{SD} = 10^{17} cm^{-3}, \it{n}','\it{N}_{SD} = 10^{17} cm^{-3}, \it{p}',...
    '\it{N}_{SD} = 0 cm^{-3}, \it{n}','\it{N}_{SD} = 0 cm^{-3}, \it{p}', 'n-type, \it{n}', 'n-type, \it{p}')
ylim([1e14, 1e20])
hold off

%% OC,ac np F x
dfplot.acnpFx(OC.intrinsic)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.acnpFx(OC.hom)
subplot(3, 1, 1)
hold off
subplot(3, 1, 2)
legend('HTL','int1','Active','int2','ETL','Intrinsic, \it{n}','Intrinsic, \it{p}','Homjunction, \it{n}','Homjunction, \it{p}')
ylim([1e14, 1e20])
hold off
subplot(3, 1, 3)
hold off

%% OC, ac np F x mobile vs static dopants
dfplot.acnpFx(OC.mobile_dope)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.acnpFx(OC.mobile_dope_i)
subplot(3, 1, 1)
legend('','','','','','','','static dopants, \it{N}_{A}', 'static dopants, \it{N}_{D}',...
    'mobile dopants, \it{N}_{A}', 'mobile dopants, \it{N}_{D}')
ylim([-0.2e17, 6e17])
hold off
subplot(3, 1, 2)
legend('','','','','','','','static dopants, \it{n}', 'static dopants, \it{p}', 'mobile dopants, \it{n}', 'mobile dopants, \it{p}')
ylim([1e15, 1e18])
hold off
subplot(3, 1, 3)
legend('ETL','interface 1','n-type','interface 2','p-type','interface 3','HTL', 'static dopants', 'mobile dopants')
ylim([-5e4, 5e4])
hold off

%% OC, ac np F x Intrinsic vs mobile vs static dopants
dfplot.rhoacnpFx2(OC.intrinsic)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoacnpFx2(OC.mobile_dope)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoacnpFx2(OC.mobile_dope_i)
subplot(3, 1, 1)
legend('ETL','interface 1','Absorber','interface 2','HTL',...
    'Mobile, Intrinsic',...
    'Static, Intrinsic',...
    'Mobile, Static Dopants',...
    'Static, Static Dopants',...
    'Mobile, Mobile Dopants',...
    'Static, Mobile Dopants')
ylim([-1e17, 6e17])
hold off
subplot(3, 1, 2)
legend('ETL','interface 1','Absorber','interface 2','HTL',...
    'Intrinsic, \it{n}','intrinsic, \it{p}',...
    'Static Dopants, \it{n}', 'Static Dopants, \it{p}',...
    'Mobile Dopants, \it{n}', 'Mobile Dopants, \it{p}')
ylim([1e15, 1e20])
hold off
subplot(3, 1, 3)
legend('ETL','interface 1','Absorber','interface 2','HTL','intrinsic','static dopants','mobile dopants')
ylim([-6e4, 6e4])
hold off

%% OC, ac np V x Intrinsic vs mobile vs static dopants
dfplot.rhoacnpFx2(OC.intrinsic)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoacnpFx2(OC.mobile_dope)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoacnpFx2(OC.mobile_dope_i)
subplot(3, 1, 1)
legend('ETL','interface 1','Absorber','interface 2','HTL',...
    'Mobile, Intrinsic',...
    'Static, Intrinsic',...
    'Mobile, Static Dopants',...
    'Static, Static Dopants',...
    'Mobile, Mobile Dopants',...
    'Static, Mobile Dopants')
ylim([-0.6e18, 1.6e18])
hold off
subplot(3, 1, 2)
legend('ETL','interface 1','Absorber','interface 2','HTL',...
    'Intrinsic, \it{n}','intrinsic, \it{p}',...
    'Static Dopants, \it{n}', 'Static Dopants, \it{p}',...
    'Mobile Dopants, \it{n}', 'Mobile Dopants, \it{p}')
ylim([1e15, 1e19])
hold off
subplot(3, 1, 3)
legend('ETL','interface 1','Absorber','interface 2','HTL','intrinsic','static dopants','mobile dopants')
ylim([-6e4, 6e4])
hold off

%% OC, npx Intrinsic vs mobile vs static dopants
dfplot.npx(OC.mobile_dope)
hold on
dfplot.npx(OC.intrinsic)
hold on
dfplot.npx(OC.mobile_dope_i)
hold off
legend('ETL','interface 1','n-type','interface 2','p-type','interface 3','HTL','static dopants, \it{n}', 'static dopants, \it{p}','intrinsic, \it{n}','intrinsic, \it{p}', 'mobile dopants, \it{n}', 'mobile dopants, \it{p}')
ylim([1e15, 1e20])

%% OC, ac np Fx Homojunciton NSD 0 1e17 1e19
dfplot.rhoacnpFx2(OC.hom_i)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoacnpFx2(OC.hom_Nion1e17_i)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoacnpFx2(OC.hom)
subplot(3, 1, 1)
legend('ETL','interface 1','n-type','interface 2','p-type','interface 3','HTL',...
    'Mobile, \it{N}_{SD} = 10^{19} cm^{-3}',...
        'Static, \it{N}_{SD} = 10^{19} cm^{-3}',...
    'Mobile, \it{N}_{SD} = 10^{17} cm^{-3}',...
        'Static, \it{N}_{SD} = 10^{17} cm^{-3}',...
        'Mobile, \it{N}_{SD} = 0 cm^{-3}',...
    'Static, \it{N}_{SD} = 0 cm^{-3}')
ylim([-0.6e18, 1.6e18])
hold off
subplot(3, 1, 2)
legend('ETL','interface 1','n-type','interface 2','p-type','interface 3','HTL',...
'\it{n}, \it{N}_{SD} = 10^{19} cm^{-3}','\it{p}, \it{N}_{SD} = 10^{19} cm^{-3}',... 
    '\it{n}, \it{N}_{SD} = 10^{17} cm^{-3}','\it{p}, \it{N}_{SD} = 10^{17} cm^{-3}',...
    '\it{n}, \it{N}_{SD} = 0 cm^{-3}', '\it{p}, \it{N}_{SD} = 0 cm^{-3}')
ylim([1e15, 1e19])
hold off
subplot(3, 1, 3)
legend('ETL','interface 1','n-type','interface 2','p-type','interface 3','HTL',...
    '\it{N}_{SD} = 10^{19} cm^{-3}','\it{N}_{SD} = 10^{17} cm^{-3}','\it{N}_{SD} = 0 cm^{-3}')
ylim([-6e4, 6e4])
hold off

%% Energy levels OC
dfplot.ELx_single(OC.intrinsic)
hold on
dfplot.ELx_single(OC.mobile_dope)
hold on
dfplot.ELx_single(OC.mobile_dope_i)
hold off

%% Energy levels rho
dfplot.rho_electronic_x(OC.intrinsic)
hold on
dfplot.rho_electronic_x(OC.mobile_dope)
hold on
dfplot.rho_electronic_x(OC.mobile_dope_i)
hold off
legend('ETL','interface 1','Absorber','interface 2','HTL',...
        'Intrinsic', 'Static Dopants', 'Mobile Dopants')

%% Energy levels SC dark
dfplot.ELx_single(soleq.intrinsic.el)
legend('ETL','interface 1','Absorber','interface 2','HTL', 'E_{fn}', 'E_{fp}', 'E_{CB}', 'E_{VB}')
ylim([-8, -1])
%%
dfplot.ELx_single(soleq.intrinsic.ion)
legend('ETL','interface 1','Absorber','interface 2','HTL', 'E_{fn}', 'E_{fp}', 'E_{CB}', 'E_{VB}')
ylim([-8, -1])
%%
dfplot.ELx_single(soleq.hom.el)
legend('ETL','interface 1','n-type','interface 2','p-type','interface 3','HTL', 'E_{fn}', 'E_{fp}', 'E_{CB}', 'E_{VB}')
ylim([-8, -1])
%%
dfplot.ELx_single(soleq.hom.ion)
legend('ETL','interface 1','n-type','interface 2','p-type','interface 3','HTL', 'E_{fn}', 'E_{fp}', 'E_{CB}', 'E_{VB}')
ylim([-8, -1])