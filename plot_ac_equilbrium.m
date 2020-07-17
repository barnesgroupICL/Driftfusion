
%% Intrinsic
dfplot.acx(soleq.intrinsic.ion)
hold on
dfplot.acx(soleq.intrinsic.el)
legend('HTL','int1','Active','int2','ETL',...
    'Intrinsic with Shottky Defects, a', 'Intrinsic with Shottky Defects, c',...
    'Intrinsic, a', 'Intrinsic, c')
hold off

%% Homojunction
dfplot.acx(soleq.hom.ion)
hold on
dfplot.acx(soleq.hom.el)
legend('HTL','int1','n-type','int2','p-type','int3','ETL',...
    'Homojunction with Shottky Defects, a','Homojunction with Shottky Defects, c',...
    'Homojunction, a', 'Homojunction, c')
hold off