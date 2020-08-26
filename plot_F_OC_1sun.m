%% Intrinsic
dfplot.Fx(OC.intrinsic_Nion1e17_i)
hold on
dfplot.Fx(OC.intrinsic_Nion1e18_i)
hold on
dfplot.Fx(OC.intrinsic_i)
hold on
dfplot.Fx(OC.intrinsic)
legend('HTL','int1','Active','int2','ETL','Intrinsic, N_{cat} = 10^{17} cm{-3}','Intrinsic, N_{cat} = 10^{18} cm{-3}',...
    'Intrinsic, N_{cat} = 10^{19} cm{-3}','Intrinsic, N_{cat} = 0 cm{-3}')
hold off
ylim([-1e5,1e5])

%% Homojunction
dfplot.Fx(OC.hom_Nion1e17_i)
hold on
dfplot.Fx(OC.hom_Nion1e18_i)
hold on
dfplot.Fx(OC.hom_i)    % Nion = 1e19 cm-3
hold on
dfplot.Fx(OC.hom)
legend('HTL','int1','n-type','int2','p-type','int3','ETL','Homojunction, N_{cat} = 10^{17} cm{-3}','Homojunction, N_{cat} = 10^{18} cm{-3}',...
   'Homojunction, N_{cat} = 10^{19} cm{-3}', 'Homojunction, N_{cat} = 0 cm{-3}')
hold off
ylim([-1e5,1e5])