function D = Ddfn(Nc, Ec, Efn, mu, T)
% Test function for Fermi Dirac diffusion coefficient
% Nc = conduction band density of states
% Ec = conduction band energy
% Ef = Fermi level
% T = temperature
% See Schubert 2015, pp. 130

kB = 8.617330350e-5;       % Boltzmann constant [eV K^-1]
kT = kB*T;
e = 1.61917e-19;         % Elementary charge in Coulombs.

for i = 1:length(Efn)
    
    f = @(E) 1./(1 + exp((E-Efn(i)+Ec)/kT));    % Fermi- Dirac function
    dfdE = @(E) exp((E-Efn(i)+Ec)/kT)./(kT*(exp((E-Efn(i)+Ec)/kT)+1).^2);
    
    g = @(E) (E/kT).^0.5;                        % DOS function based on 3D semiconductor
    h = @(E) g(E).*f(E);
    k = @(E) g(E).*dfdE(E);
    
    n(i) = real(((2*Nc)/(kT*pi^0.5))*integral(f, 0, 16));
    dndE(i) = ((2*Nc)/(kT*pi^0.5))*integral(dfdE, 0, 16);
    
end

D = mu*(n./dndE);

figure(1)
semilogy(Efn, n)
xlabel('Fermi energy [eV]')
ylabel('Carrier density [cm-3]')
hold on

figure(2)
semilogy(Efn, dndE)
xlabel('Fermi energy [eV]')
ylabel('dndEfn [cm-3/eV]')
hold on

figure(3)
semilogy(Efn, D)
xlabel('Fermi energy [eV]')
ylabel('Diffusion coefficient [cm2s-1]')
hold on

figure(4)
loglog(n, D)
xlabel('Carrier density [cm-3]')
ylabel('Diffusion coefficient [cm2s-1]')
hold on
end