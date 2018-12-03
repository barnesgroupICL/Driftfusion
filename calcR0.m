function [JV_ana, R0, k_rad, Voc, G0] = calcR0(Eg, N0, EgArr, Jsc_vs_Eg, par)
% calculates 
% use calcJsc first to get Jsc_vs_Eg & egArr
%% Input arguments
% PAR - a parameters class containing 

set(0,'DefaultLineLinewidth',1);
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
set(0,'DefaultAxesXcolor', [0, 0, 0]);
set(0,'DefaultAxesYcolor', [0, 0, 0]);
set(0,'DefaultAxesZcolor', [0, 0, 0]);
set(0,'DefaultTextColor', [0, 0, 0]);

figson = 1;

%% Physical constants
h = 4.135667662e-15;        % Plancks constant [eVs-1]
c = 29979245800;            % Speed of light [cms-1]
kB = 8.617330350e-5;        % Boltzmann constant [eV K^-1]
q = 1.60217662e-19;         % Elementary charge [C]

% E_dis = 1.6:0.01:1600;
% r0_dis = ((2*pi)/(h^3*c^2))*((E_dis.^2)./(exp((E_dis/(par.kB*par.T))-1)));
% 
% R01 = trapz(E_dis, r0_dis);
% 
% figure(499)   
% semilogy(E_dis, r0_dis)
% xlabel('Energy [eV]')
% ylabel('r0(E)')

for i =1:length(par.Eg)

% Find maximum Jsc based on step function absorption and 100% EQE
p = find(EgArr <= par.Eg(i));
p = p(end);

Jsc(i) = Jsc_vs_Eg(p);

fun = @(E) ((2*pi)/(h^3*c^2))*((E.^2)./(exp((E/(par.kB*par.T))-1)));

J0_fd(i) = integral(fun, Eg(i), Inf);    % Spectral bb flux density- factor of 2 for back reflector
R0(i) = J0_fd./d(i);

J0(i) = J0_fd*q;

k_rad(i) = R0/(par.ni(i)^2);      % cm3s-1

V = -0.4:0.001:1.2;

J = (-Jsc + J0*exp((V/(par.kB*par.T))-1));  %mA cm-2

J_dk = (J0*exp((V/(par.kB*par.T))-1));

Voc = (par.kB*par.T)*(log(Jsc/J0)+1);

G0 = Jsc./(t*q);           % Generation rate cm-3s-1

JV_ana.J = J;
JV_ana.V = V;

if figson == 1

figure(500)
plot(V, J, V, J_dk)
xlabel('V [V]')
ylabel('J[mAcm-2]')
ylim([-50, 30])
xlim([-inf, inf])
grid off

end

end