%% npx
dfplot.npx(JV.intrinsic.ill.f, 0.9939e3)
hold on
dfplot.npx(JV.intrinsic_i.ill.f, 1.018e3)
ylim([1e15, 1e18])
xlim([0, 600])
legend('','','','','','n, NSD = 0 cm-3', 'p, NSD = 0 cm-3',...
    'n, NSD = 1019 cm-3', 'p, NSD = 1019 cm-3')

%% acx
dfplot.acx(JV.intrinsic_i.ill.f, 1.018e3)
xlim([0, 600])