function [JV_ana, R0, k_rad, Voc, G0] = calcR0(Eg, t, EgArr, Jsc_vs_Eg)

params = pc;

% use calJsc first to get Jsc_vs_Eg & egArr
% t = thickness in cm

set(0,'DefaultLineLinewidth',1);
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
set(0,'DefaultAxesXcolor', [0, 0, 0]);
set(0,'DefaultAxesYcolor', [0, 0, 0]);
set(0,'DefaultAxesZcolor', [0, 0, 0]);
set(0,'DefaultTextColor', [0, 0, 0]);

figson = 1;

h = 4.135667662e-15;        % eVs-1
c = 29979245800;            % cms-1
Vth = params.kB*params.T;               % eV
q = params.e;         % C

%t = 400e-7;%params.xmax;                % cm
%E = 1.6:0.01:100;

% E_dis = 1.6:0.01:1600;
% r0_dis = ((2*pi)/(h^3*c^2))*((E_dis.^2)./(exp((E_dis/Vth)-1)));
% 
% R01 = trapz(E_dis, r0_dis);
% 
% figure(499)   
% semilogy(E_dis, r0_dis)
% xlabel('Energy [eV]')
% ylabel('r0(E)')

% Find maximum Jsc based on step function absorption and 100% EQE
p = find(EgArr <= Eg);
p = p(end);

Jsc = Jsc_vs_Eg(p);

fun = @(E) ((2*pi)/(h^3*c^2))*((E.^2)./(exp((E/Vth)-1)));

J0_fd = integral(fun, Eg, Inf);    % Spectral bb flux density- factor of 2 for back reflector
R0 = J0_fd./t;

J0 = J0_fd*q*1e3;

ni = params.ni;             % cm-3

k_rad = R0/(ni(2)^2);      % cm3s-1

V = -0.4:0.001:1.2;

J = (-Jsc + J0*exp((V/Vth)-1));  %mA cm-2

J_dk = (J0*exp((V/Vth)-1));

Voc = Vth*(log(Jsc/J0)+1);

G0 = (Jsc*1e-3)/(t*q);           % Generation rate cm-3s-1

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