function blakemore_boltz(Nc, eta, gamma)
% ETA is an array containing reduced electron Fermi energies (Efn -
% Ecb)/kB*T

n_blake = Nc./(exp(-eta) + gamma);
n_boltz = Nc./exp(-eta);

eta_blake_check = log(n_blake./(Nc - gamma*n_blake));
eta_boltz_check = log(n_boltz./Nc);

g_blake = (Nc./(Nc - gamma*n_blake));
g_boltz = ones(1, length(eta));

figure(1)
semilogy(eta, n_blake, eta, n_boltz)
xlabel('(Efn - Ecb)/kBT')
ylabel('n [cm-3]')
legend('Blakemore', 'Boltz')

figure(2)
semilogy(eta, g_blake, eta, g_boltz)
xlabel('(Efn - Ecb)/kBT')
ylabel('g(epsilon)')
legend('Blakemore', 'Boltz')

figure(3)
plot(eta, eta_blake_check, eta, eta_boltz_check)
xlabel('(Efn - Ecb)/kBT')
ylabel('Check')
legend('Blakemore', 'Boltz')

end