%% 
dfplot.JV(JV.intrinsic, 2)
figure(4)
hold on
dfplot.JV(JV.intrinsic_i, 2)
figure(4)
hold on
dfplot.JV(JV.hom, 2)
figure(4)
hold on
dfplot.JV(JV.hom_i, 2)
hold on
dfplot.JV(JV.ntype_i, 2)
figure(4)
hold off
legend('', 'Intrinsic', '', 'Intrinsic with Shottky Defects',...
    '','Homojunction', '','Homojunction with Shottky Defects', '','n-type')
ylim([-5e-3, 25e-3])

