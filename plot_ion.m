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

%% Equilbrium EL
dfplot.ELx_single(soleq.intrinsic.ion)
hold on
dfplot.ELx_single(soleq.ntype.ion)
hold on
dfplot.ELx_single(soleq.hom.ion)
hold on
dfplot.ELx_single(soleq.mobile_dope.ion)
legend('HTL','int1','Active','int2','ETL','intrinsic','ntype','homojunction','mobile dopants')
hold off

%% rhox Fx Vx
dfplot.rhoxFxVx(OC.intrinsic_i)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC.ntype_i)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC.hom_i)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC.mobile_dope_i)
subplot(3, 1, 1)
legend('HTL','int1','Active','int2','ETL','intrinsic','ntype','homojunction','mobile dopants')
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on