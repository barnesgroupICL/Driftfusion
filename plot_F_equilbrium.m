
%% Intrinsic
dfplot.Fx(soleq.intrinsic.ion)
hold on
dfplot.Fx(soleq.intrinsic.el)
legend('HTL','int1','Active','int2','ETL','Intrinsic with Shottky Defects','Intrinsic')
hold off

%% Homojunction
dfplot.Fx(soleq.hom.ion)
hold on
dfplot.Fx(soleq.hom.el)
legend('HTL','int1','n-type','int2','p-type','int3','ETL','Homojunction with Shottky Defects','Homojunction')
hold off