function [JV_ana, r0, Voc, g0] = calcR0_pn_junction(par)
% calculates
% use calcJsc first to get Jsc_vs_Eg & egArr
%% Input arguments
% PAR - a parameters class containing
figson = 1;

%% Physical constants
h = 4.135667662e-15;        % Plancks constant [eVs-1]
c = 29979245800;            % Speed of light [cms-1]
kB = 8.617330350e-5;        % Boltzmann constant [eV K^-1]
q = par.e;%1.60217662e-19;         % Elementary charge [C]

% E_dis = 1.6:0.01:1600;
% r0_dis = ((2*pi)/(h^3*c^2))*((E_dis.^2)./(exp((E_dis/(par.kB*par.T))-1)));
%
% R01 = trapz(E_dis, r0_dis);
%
% figure(499)
% semilogy(E_dis, r0_dis)
% xlabel('Energy [eV]')
% ylabel('r0(E)')
[EgArr,Jsc_vs_Eg] = calcJsc(0);

%for i =1:length(par.Eg)
    
    % Find maximum Jsc based on step function absorption and 100% EQE
    p = find(EgArr <= par.Eg(1));
    p = p(end);
    
    Jsc = Jsc_vs_Eg(p);
    
    fun = @(E) ((2*pi)/(h^3*c^2))*((E.^2)./(exp((E/(par.kB*par.T))-1)));
    
    j0 = integral(fun, par.Eg(1), Inf);    % Spectral bb flux density- factor of 2 for back reflector
    
    r0 = j0./par.dcum(end);
    
    taun = par.taun(1);
    taup = par.taup(end);
    
    Dn = par.mue(1)*par.kB*par.T;
    Dp = par.muh(1)*par.kB*par.T;
    Ln = (taun*Dn)^0.5;
    Lp = (taup*Dp)^0.5;
    J0 = (q*Dp*par.p0(end)/Lp) + (q*Dn*par.n0(1)/Ln);
    
    %k_rad(i) = r0(i)/(par.ni(i)^2);      % cm3s-1
    
    V = 0:1e-3:1.2;
    
    J = -Jsc + J0.*(exp(V./(par.kB*par.T))-1);  %A cm-2
    
    J_dk = J0.*(exp(V./(par.kB*par.T))-1);
    
    Voc = (par.kB*par.T).*(log(Jsc./J0)+1);
    
              % Generation rate cm-3s-1
    
    JV_ana.J = J;
    JV_ana.J_dk = J_dk;
    JV_ana.V = V;
    
    if figson == 1
        figure(4)
        plot(V, J_dk, V, J)
        xlabel('Applied Voltage [V]')
        ylabel('Current density [Acm-2]')
        ylim([-50e-3, 30e-3])
        %xlim([0, V(end)])
    end
    
%end
    g0 = Jsc./(par.dcum(end).*q); 

end