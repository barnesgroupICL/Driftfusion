function D = Ddfn_numeric(Nc, Ec, Efn, mu, T)
% Test function for Fermi Dirac diffusion coefficient 
% Nc = conduction band density of states
% Ec = conduction band energy
% Ef = Fermi level
% T = temperature
% See Schubert 2015, pp. 130

kB = 8.617330350e-5;       % Boltzmann constant [eV K^-1]    
kT = F.kB*T;
e = 1.61917e-19;         % Elementary charge in Coulombs.

for i = 1:length(Efn)
    
    E = Ec:0.01:(10+Ec);
    f = 1./(1 + exp((E-Efn(i)+Ec)/kT));    % Fermi- Dirac function
    dfdE = ((exp((E-Efn(i)+Ec)/kT))./(kT*((exp((E-Efn(i)+Ec)/kT)+1).^2)));

    g = (E/kT).^0.5;                        % DOS function based on 3D semiconductor
    h = g.*f;
    k = g.*dfdE;

    n(i) = real(((2*Nc)/(kT*pi^0.5))*trapz(E, h));
    dndE(i) = real(((2*Nc)/(kT*pi^0.5))*trapz(E, k));
end

D = mu*(n./dndE);

% Plotting defaults
set(0,'DefaultLineLinewidth',1);
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
set(0,'DefaultAxesXcolor', [0, 0, 0]);
set(0,'DefaultAxesYcolor', [0, 0, 0]);
set(0,'DefaultAxesZcolor', [0, 0, 0]);
set(0,'DefaultTextColor', [0, 0, 0]);

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
%ylim([0, 0.05])
end