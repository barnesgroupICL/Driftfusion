%% Without ions
dfplot.JV(JV.standard_epp5, 2)
figure(4) 
hold on
dfplot.JV(JV.standard_epp25, 2)
figure(4) 
hold on
dfplot.JV(JV.standard_epp50, 2)
figure(4) 
hold on
dfplot.JV(JV.standard_i, 2)
figure(4) 
hold off
xlim([0,1.2])
ylim([-30e-3, 10e-3])
legend('epp = 5, F', 'epp = 5, R', 'epp = 25, F', 'epp = 25, R',...
    'epp = 50 F', 'epp = 50 R', 'NSD = 1e19cm-3, F', 'NSD = 1e19cm-3, F')

%% With ions
dfplot.JV(JV.pin_srh_i, 2)
figure(4) 
hold on
dfplot.JV(JV.hom_srh_i, 2)
figure(4) 
hold on
dfplot.JV(JV.ntype_srh_i, 2)
figure(4) 
hold on
dfplot.JV(JV.mobile_dope_i, 2)
figure(4) 
hold on