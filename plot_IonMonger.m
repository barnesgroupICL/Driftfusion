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
    
%% Profiles
p_array = [1, 21, 41, 61, 81];

sol_im = sol_light_noIR_1Vs_im;
x_im = [sol_im.vectors.xE; sol_im.vectors.x; sol_im.vectors.xH]+100;
n_im = [sol_im.dstrbns.nE, sol_im.dstrbns.n, zeros(length(sol_im.time), size(sol_im.dstrbns.pH, 2))]*1e-6;
p_im = [zeros(length(sol_im.time), size(sol_im.dstrbns.nE, 2)), sol_im.dstrbns.p, sol_im.dstrbns.pH]*1e-6;
c_im = [zeros(length(sol_im.time), size(sol_im.dstrbns.nE, 2)),  sol_im.dstrbns.P, zeros(length(sol_im.time), size(sol_im.dstrbns.pH, 2))]*1e-6;    
V_im = [sol_im.dstrbns.phiE, sol_im.dstrbns.phi, sol_im.dstrbns.phiH]; 

figure(13); hold on
figure(14); hold on
for i =1:length(p_array)
    % electron profile
    figure(12)
    plot(x_im, V_im(i, :), 'k--')
    
    figure(13)
    plot(x_im, n_im(i, :), 'k--')
    plot(x_im, p_im(i, :), 'k--')
    
    figure(14)
    plot(x_im, c_im(i, :), 'k--')
end
    