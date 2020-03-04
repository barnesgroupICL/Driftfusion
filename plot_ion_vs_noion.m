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
legend('HTL','int1','Active','int2','ETL','intrinsic NSD=0cm-3','intrinsic NSD=0cm-3')
subplot(2, 1, 2)
hold off

%%
dfplot.acnpFx(OC.mobile_dope)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.acnpFx(OC.mobile_dope_i)
subplot(3, 1, 1)
legend('','','','','','','','static dopants, NA', 'static dopants, ND',...
    'mobile dopants, NA', 'mobile dopants, ND')
xlim([0, 600])
ylim([-0.2e17, 6e17])
hold off
subplot(3, 1, 2)
legend('','','','','','','','static dopants, n', 'static dopants, p', 'mobile dopants, n', 'mobile dopants, p')
xlim([0, 600])
ylim([1e15, 1e18])
hold off
subplot(3, 1, 3)
legend('ETL','interface 1','n-type','interface 2','p-type','interface 3','HTL', 'static dopants', 'mobile dopants')
xlim([0, 600])
ylim([-5e4, 5e4])
hold off