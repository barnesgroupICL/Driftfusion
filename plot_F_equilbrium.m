
%% Intrinsic
dfplot.Fx(soleq.intrinsic_Nion1e17.ion)
hold on
dfplot.Fx(soleq.intrinsic_Nion1e18.ion)
hold on
dfplot.Fx(soleq.intrinsic.ion)
hold on
dfplot.Fx(soleq.intrinsic.el)
legend('HTL','int1','Active','int2','ETL','Intrinsic, N_{cat} = 10^{17} cm{-3}','Intrinsic, N_{cat} = 10^{18} cm{-3}',...
    'Intrinsic, N_{cat} = 10^{19} cm{-3}','Intrinsic, N_{cat} = 0 cm{-3}')
hold off
ylim([-0.5e5,2e5])

%% Homojunction
dfplot.Fx(soleq.hom_Nion1e17.ion)
hold on
dfplot.Fx(soleq.hom_Nion1e18.ion)
hold on
dfplot.Fx(soleq.hom.ion)    % Nion = 1e19 cm-3
hold on
dfplot.Fx(soleq.hom.el)
legend('HTL','int1','n-type','int2','p-type','int3','ETL','Homojunction, N_{cat} = 10^{17} cm{-3}','Homojunction, N_{cat} = 10^{18} cm{-3}',...
    'Homojunction N_{cat} = 10^{19} cm{-3}', 'Homojunction N_{cat} = 0 cm{-3}')
hold off
ylim([-0.5e5,2e5])
