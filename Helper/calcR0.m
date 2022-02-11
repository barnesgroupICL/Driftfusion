function [JV_ana, r0, k_rad, Voc, g0] = calcR0(EgArr, Jsc_vs_Eg, par)
% calculates
% use calcJsc first to get Jsc_vs_Eg & egArr
% Input arguments
% PAR - a parameters class containing
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Start code
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

    if isnan(par.Eg(i))
        continue
    end

    % Find maximum Jsc based on step function absorption and 100% EQE
    p = find(EgArr <= par.Eg(i));
    p = p(end);

    Jsc(i) = Jsc_vs_Eg(p);

    fun = @(E) ((2*pi)/(h^3*c^2))*((E.^2)./(exp((E/(par.kB*par.T))-1)));

    j0(i) = integral(fun, par.Eg(i), Inf);    % Spectral bb flux density- factor of 2 for back reflector
    r0(i) = j0(i)./par.d(i);

    J0(i) = j0(i)*q;

    k_rad(i) = r0(i)/(par.ni(i)^2);      % cm3s-1

    V = -0.4:1e-2:3.0;

    J(i,:) = -Jsc(i) + J0(i).*(exp(V./(par.kB*par.T))-1);  %A cm-2

    J_dk(i,:) = J0(i).*(exp(V/(par.kB*par.T))-1);

    Voc(i) = (par.kB*par.T).*(log(Jsc(i)./J0(i))+1);

    g0(i) = Jsc(i)./(par.d(i).*q);           % Generation rate cm-3s-1

    JV_ana.J(i,:) = J(i,:);
    JV_ana.J_dk(i,:) = J_dk(i,:);
    JV_ana.V = V;

    if figson == 1
        figure(500)
        plot(V, J(i,:), V, J_dk(i,:))
        xlabel('V [V]')
        ylabel('J[mAcm-2]')
        ylim([-50, 30])
        xlim([0, V(end)])
        hold on
    end
end
    figure(500)
    hold off
end
