%% Interfacial recombination
figure(91)
hold on
plotJV_im(sol_light_50mVs);
plotJV_im(sol_light_100mVs);
plotJV_im(sol_light_200mVs);
ylim([-10e-3,25e-3])
xlim([-1.2, 0])
legend('DF- 50 mVs-1', 'DF- 100 mVs-1', 'DF- 200 mVs-1',...
        'IM- 50 mVs-1', 'IM- 100 mVs-1', 'IM- 200 mVs-1')

%% Bulk recombination
figure(91)
hold on
plotJV_im(sol_light_noIR_0p1Vs_im);
plotJV_im(sol_light_noIR_1Vs_im);
plotJV_im(sol_light_noIR_10Vs_im);
ylim([-10e-3,25e-3])
xlim([-1.2, 0])
legend('DF- 0.1 Vs-1', 'DF- 1 Vs-1', 'DF- 10 Vs-1',...
        'IM- 0.1 Vs-1', 'IM- 1 Vs Vs-1', 'IM- 10 Vs-1')
    
    