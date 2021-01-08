%% JV plots
figure(91)
hold on
% plotJV(sol_dark_50mVs);
% plotJV(sol_dark_100mVs);
% plotJV(sol_dark_200mVs);
plotJV_im(im_sol_light_50mVs);
plotJV_im(im_sol_light_100mVs);
plotJV_im(im_sol_light_200mVs);
ylim([-10e-3,25e-3])
xlim([0, 1.1])
legend('DF 50 mVs-1', 'DF 100 mVs-1', 'DF 200 mVs-1',...
    'IM 50 mVs-1', 'IM 100 mVs-1', 'IM 200 mVs-1') 