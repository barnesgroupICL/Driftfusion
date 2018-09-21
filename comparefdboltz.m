function comparefdboltz
% Plotting defaults
set(0,'DefaultLineLinewidth',1);
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
set(0,'DefaultAxesXcolor', [0, 0, 0]);
set(0,'DefaultAxesYcolor', [0, 0, 0]);
set(0,'DefaultAxesZcolor', [0, 0, 0]);
set(0,'DefaultTextColor', [0, 0, 0]);

Ef = -0.3:0.01:0.1;

Nc = 1e20;
Ec = 0;
T = 300;

for i = 1:length(Ef)
    
    n_fd(i) = F.fdn(Nc, Ec, Ef(i), T);
    
    n_boltz(i) = F.boltzn(Nc, Ec, Ef(i), T);
    
end

semilogy((Ef-Ec), n_fd,  (Ef-Ec), n_boltz)
xlabel('Fermi energy [eV]')
ylabel('Carrier density [cm-3]')
legend('Fermi-Dirac', 'Boltzmann')

end