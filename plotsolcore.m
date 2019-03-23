function plotsolcore(solcoresol)

% Plotting defaults
set(0,'DefaultLineLinewidth',1);
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
set(0,'DefaultAxesXcolor', [0, 0, 0]);
set(0,'DefaultAxesYcolor', [0, 0, 0]);
set(0,'DefaultAxesZcolor', [0, 0, 0]);
set(0,'DefaultTextColor', [0, 0, 0]);

s = solcoresol;
xnm = s.x*1e9;

figure(1)
subplot(3, 1, 1)
hold on
plot(xnm, s.Ec, xnm, s.Ev, xnm, s.Efn, '--', xnm, s.Efp, '--')
xlabel('Position')
ylabel('Energy [eV]')
xlim([0, xnm(end)])

subplot(3, 1, 2)
hold on
semilogy(xnm, s.n*1e-6, xnm, s.p*1e-6)
xlabel('Position')
ylabel('Carrier density [cm-3]')
xlim([0, xnm(end)])

figure(251)
plot(xnm, s.EA, xnm, s.IP)
legend('EA', 'IP')
xlabel('Position')
ylabel('Band energies [eV]')
xlim([0, xnm(end)])

figure(252)
plot(xnm, s.V)
xlabel('Position')
ylabel('El Potential [V]')
xlim([0, xnm(end)])

end