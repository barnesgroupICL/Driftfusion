%% Without ions
dfplot.JV(JV.intrinsic, 2)
figure(4)
hold on
dfplot.JV(JV.ntype, 2)
figure(4)
hold on
dfplot.JV(JV.hom, 2)
figure(4)
hold on
dfplot.JV(JV.mobile_dope, 2)
figure(4)
hold off
legend('intrinsic', 'homojunciton', 'n-type', 'mobile dopants')
ylim([-5e-3, 25e-3])

%% With ions
dfplot.JV(JV.intrinsic_i, 2)
figure(4)
hold on
dfplot.JV(JV.hom_i, 2)
figure(4)
hold on
dfplot.JV(JV.ntype_i, 2)
figure(4)
hold on
dfplot.JV(JV.mobile_dope_i, 2)
figure(4)
hold off
legend('intrinsic', 'homojunciton', 'n-type', 'mobile dopants')
ylim([-5e-3, 25e-3])