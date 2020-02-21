%% JVs
dfplot.JV(JV.symmetric_srh_epp5, 2)
figure(4) 
hold on
dfplot.JV(JV.symmetric_srh_epp25, 2)
figure(4) 
hold on
dfplot.JV(JV.symmetric_srh_epp50, 2)
figure(4) 
hold on
dfplot.JV(JV.symmetric_srh_epp25_i, 2)
figure(4) 
hold off
xlim([0,1.2])
ylim([-30e-3, 10e-3])
legend('epp = 5, F', 'epp = 5, R', 'epp = 25, F', 'epp = 25, R',...
    'epp = 50, F', 'epp = 50, R', 'NSD = 1e19cm-3, F', 'NSD = 1e19cm-3, R')

%% OC vs time
dfplot.Voct(OC_Newton.symmetric_srh_epp5);
hold on
dfplot.Voct(OC_Newton.symmetric_srh_epp25);
hold on
dfplot.Voct(OC_Newton.symmetric_srh_epp50);
hold on
dfplot.Voct(OC_Newton.symmetric_srh_epp25_i);
hold off

%% OC EL
dfplot.ELx(OC_Newton.symmetric_srh_epp5);
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.ELx(OC_Newton.symmetric_srh_epp25);
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.ELx(OC_Newton.symmetric_srh_epp50);
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.ELx(OC_Newton.symmetric_srh_epp25_i);
subplot(3, 1, 1)
hold off
subplot(3, 1, 2)
hold off
subplot(3, 1, 3)
hold off

%% rhox Fx Vx
dfplot.rhoxFxVx(OC_Newton.symmetric_srh_epp25_i);
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC_Newton.symmetric_srh_epp5);
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC_Newton.symmetric_srh_epp25);
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.rhoxFxVx(OC_Newton.symmetric_srh_epp50);
subplot(3, 1, 1)
legend('HTL','int1','Active','int2','ETL',...
    'NSD = 1e19 cm-3, epp = 25',...
    'NSD = 0 cm-3, epp = 5',...
    'NSD = 0 cm-3, epp = 25',...
    'NSD = 0 cm-3, epp = 50')
hold off
subplot(3, 1, 2)
hold off
subplot(3, 1, 3)
hold off

%% rhox
dfplot.rhox(OC_Newton.symmetric_srh_epp25_i);
hold on
dfplot.rhox(OC_Newton.symmetric_srh_epp5);
hold on
dfplot.rhox(OC_Newton.symmetric_srh_epp25);
hold on
dfplot.rhox(OC_Newton.symmetric_srh_epp50);
hold off
legend('','','','','','ion, epp = 25', 'epp = 5', 'epp = 25', 'epp = 50')

%% Vx
dfplot.Vx(OC_Newton.symmetric_srh_epp25_i);
hold on
dfplot.Vx(OC_Newton.symmetric_srh_epp5);
hold on
dfplot.Vx(OC_Newton.symmetric_srh_epp25);
hold on
dfplot.Vx(OC_Newton.symmetric_srh_epp50);
hold off
legend('','','','','','ion, epp = 25', 'epp = 5', 'epp = 25', 'epp = 50')
