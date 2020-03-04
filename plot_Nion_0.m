%% JVs
dfplot.JV(JV.intrinsic, 2)
figure(4)
ylim([-10e-3, 30e-3])
hold on
dfplot.JV(JV.ntype, 2)
figure(4)
hold on
dfplot.JV(JV.hom, 2)
figure(4)
hold on
dfplot.JV(JV.mobile_dope, 2)
figure(4)
hold on

%% Equilbrium EL
dfplot.ELx(soleq.intrinsic.el)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.ELx(soleq.ntype.el)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.ELx(soleq.hom.el)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.ELx(soleq.mobile_dope.el)
subplot(3, 1, 1)
legend('HTL','int1','Active','int2','ETL','intrinsic','ntype','homojunction','mobile dopants')
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on

%% rhox Fx Vx
dfplot.rhoxFxVx(OC.intrinsic)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC.ntype)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC.hom)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC.mobile_dopants)
subplot(3, 1, 1)
legend('HTL','int1','Active','int2','ETL','intrinsic','ntype','homojunction','mobile dopants')
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on