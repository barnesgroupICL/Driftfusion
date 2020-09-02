%% Additional plots for comparison with ASA
% Note: Data stored locally- not available in the repository

%% Absolute error fig 10b
semilogy(DF_Vapp, abs(dk_error_percent_m1a),...
    DF_Vapp, abs(dk_error_percent_m1b),...
    DF_Vapp, abs(dk_error_percent_m2a),...
    DF_Vapp, abs(dk_error_percent_m2b))
xlabel('Applied voltage [V]')
ylabel('Absolute difference [%]')
xlim([0, 1.4])
%ylim([0, 20])
legend('PS 1a', 'PS 1b', 'PS 2a', 'PS 2b')

%% Equilibrium energy level diagrams
dfplot.ELx(soleq.no_ion)
