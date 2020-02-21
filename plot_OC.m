%% rhox Fx Vx
dfplot.rhoxFxVx(OC.intrinsic_i)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC.intrinsic_epp5)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC.intrinsic)
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC.intrinsic_epp50)
subplot(3, 1, 1)
legend('HTL','int1','Active','int2','ETL',...
    'NSD = 1e19 cm-3, epp = 25',...
    'NSD = 0 cm-3, epp = 5',...
    'NSD = 0 cm-3, epp = 25',...
    'NSD = 0 cm-3, epp = 50')
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on

%% OC vs time
dfplot.Voct(OC.intrinsic_epp5);
hold on
dfplot.Voct(OC.intrinsic_epp25);
hold on
dfplot.Voct(OC.intrinsic_epp50);
hold on
dfplot.Voct(OC.intrinsic_i);
hold off