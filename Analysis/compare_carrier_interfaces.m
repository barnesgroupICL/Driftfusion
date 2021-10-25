function [n_ana, p_ana, jn_ana, jp_ana] = compare_carrier_interfaces(sol, tarr, plot_switch)
% Function to compare the interfacial carrier densities between simulation
% and analytical solution
%
%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Version X2
% Version history
% v5 Calculated on the whole mesh, uses jn instead of js
% substitutes r for r=(js-j)/x
% v6 Tests expressions for ns, ps for majority carrier sides of the
% v7 as with v6 but calculated on the half mesh
% v8 uses left hand boundary values only and adapts expression for n(d) =
% ns and p(d) = ps dependent on sign of alpha and beta
% v8.2 includes a switch to plot ns ps calculated from n(x), p(x), jn(x),
% jp(x), jns and jps
% vX As with 8.2 but with more plot options
% vX2 Generalised for any number of interfaces and recombination zone
% TARR = array of time points
plotswitch = 1;     % 1 = plot n(x), p(x), 2 = plot ns, ps, 3 = 4x panels
calc_option = 2;    %1 = all terms of expression, 2 = first term only

par = sol.par;
T = par.T;
kB = par.kB;
dev = par.dev_sub;
u = sol.u;
t = sol.t;
x = sol.x;
x_sub = par.x_sub;
pcum1 = par.pcum + 1;   % Cumulative layer points array

xprime = dev.xprime;
alpha0 = dev.alpha0;
beta0 = dev.beta0;
mu_n = dev.mu_n;
mu_p = dev.mu_p;
dint = dev.dint;
Vt = kB*T;
kB = par.kB;
T = par.T;

dVdx = zeros(length(t), length(x_sub));
n = zeros(length(t), length(x_sub));
p = zeros(length(t), length(x_sub));
n_ana = NaN(length(t), length(x_sub));
p_ana = NaN(length(t), length(x_sub));
jn_ana = NaN(length(t), length(x_sub));
jp_ana = NaN(length(t), length(x_sub));

%% Avoid PDEVAL for faster calculation
% Obtain variables and gradients on sub-interval mesh
for i = 1:length(t)
    n(i,:) = 0.5*(u(i, 2:end, 2) + u(i, 1:end-1, 2));
    p(i,:) = 0.5*(u(i, 2:end, 3) + u(i, 1:end-1, 3));
    dVdx(i,:) = (u(i, 2:end, 1) - u(i, 1:end-1, 1))./(x(2:end) - x(1:end-1));
end
alpha = par.q*dVdx./(kB*T) + alpha0;
beta = par.q*-dVdx./(kB*T) + beta0;

% Fluxes - approximate as on half mesh
[~, j, ~] = dfana.calcJ(sol);
jn = j.n;
jp = j.p;

% recombination
r = dfana.calcr(sol, "sub");

%% Get indexes of interfaces
%int_index = find(contains(par.layer_type, 'interface'));   % only
%compatible from 2016 onwards
int_logical = strcmp(par.layer_type, 'interface');
loc = find(int_logical); % interface layer locations

ns = zeros(length(t), length(loc));     % Store time array of each ns in new column
ps = zeros(length(t), length(loc));     % Store time array of each ps in new column
jns = zeros(length(t), length(loc));     % Store time array of each ns in new column
jps = zeros(length(t), length(loc));     % Store time array of each ps in new column
r_vsr = zeros(length(t), length(x_sub));

for m = 1:length(loc)
    p_L = pcum1(loc(m)-1);
    p_R = pcum1(loc(m));
    % Check location of interfacial surface carrier densities
    ns(:, m) = n(:, p_L);
    jns(:, m) = jn(:, p_L);     % jns and jps always from LH boundary
    ps(:, m) = p(:, p_L);
    jps(:, m) = jp(:, p_L);      % jns and jps always from LH boundary
    r_vsr(:, p_L:p_R) = cumtrapz(x_sub(p_L:p_R), r.vsr(:, p_L:p_R), 2);
    
    for kk = 1:length(tarr)
        k = find(sol.t <= tarr(kk));
        k = k(end);                     % k is the time point
        for i = 1:length(x_sub)
            if x_sub(i) >= par.dcum(loc(m)-1) && x_sub(i) <= par.dcum(loc(m))
                jn_ana(k,i) = jns(k,m) - r_vsr(k,i);
                jp_ana(k,i) = jps(k,m) - r_vsr(k,i);
                %% Carrier densities
                switch calc_option
                    case 1  % Full expressions
                        if alpha(k,i) <= 0
                            X = exp(alpha(k,i)*xprime(i));
                            n_ana(k,i) = ns(k,m)*X + (jns(k,m)/(alpha(k,i)*mu_n(i)*Vt))*(1 - X)...
                                - ((jns(k,m)-jn(k,i))/(alpha(k,i)^2*mu_n(i)*Vt*xprime(i)))*(1 - X + alpha(k,i)*xprime(i));
                            
                        elseif alpha(k,i) > 0
                            Xstar = exp(alpha(k,i)*(xprime(i)-dint(k)));
                            n_ana(k,i) = ns(k,m)*Xstar + (jns(k,m)/(alpha(k,i)*mu_n(i)*Vt))*(1 - Xstar)...
                                - ((jns(k,m)-jn(k,i))/(alpha(k,i)^2*mu_n(i)*Vt*xprime(i)))*(1 + alpha(k,i)*xprime(i) - (1 + alpha(k,i)*dint(k))*Xstar);
                        end
                        
                        if beta(k,i) <= 0
                            Y = exp(beta(k,i)*xprime(i));
                            p_ana(k,i) = ps(k,m)*Y + (jps(k,m)/(beta(k,i)*mu_p(i)*Vt))*(1 - Y)...
                                - ((jps(k,m)-jp(k,i))/(beta(k,i)^2*mu_p(i)*Vt*xprime(i)))*(1 - Y + beta(k,i)*xprime(i));
                            
                        elseif beta(k,i) > 0
                            Ystar = exp(beta(k,i)*(xprime(i)-dint(k)));
                            p_ana(k,i) = ps(k,m)*Ystar + (jps(k,m)/(beta(k,i)*mu_p(i)*Vt))*(1 - Ystar)...
                                - ((jps(k,m)-jp(k,i))/(beta(k,i)^2*mu_p(i)*Vt*xprime(i)))*(1 + beta(k,i)*xprime(i) - (1 + beta(k,i)*dint(k))*Ystar);
                        end
                        %% Fluxes
                    case 2  % First terms only
                        if alpha(k,i) <0
                            X = exp(alpha(k,i).*xprime(i));
                            n_ana(k,i) = ns(k,m)*X;
                            
                        elseif alpha(k,i) > 0
                            Xstar = exp(alpha(k,i)*(xprime(i)-dint(k)));
                            n_ana(k,i) = ns(k,m)*Xstar;
                        end
                        
                        if beta(k,i) <0
                            Y = exp(beta(k,i)*xprime(i));
                            p_ana(k,i) = ps(k,m)*Y;
                            
                        elseif beta(k,i) >= 0
                            Ystar = exp(beta(k,i)*(xprime(i)-dint(k)));
                            p_ana(k,i) = ps(k,m)*Ystar;
                        end
                end
            end
        end
    end
end

if plot_switch
    %% Plotting
    for m = 1:length(loc)
        if plotswitch == 1
            %% Plot carrier densities
            figure(600)
            subplot(1,length(loc),m)
            dfplot.colourblocks(sol, [1e-20,1e20]);
            set(gca, 'YScale','log');
            for i = 1:length(tarr)
                k = find(sol.t <= tarr(i));
                k = k(end);
                
                semilogy(x_sub*1e7, n(k, :), x_sub*1e7, p(k, :),...
                    x_sub*1e7, n_ana(k,:), 'k-.', x_sub*1e7, p_ana(k,:), 'k--')
                hold on
            end
            hold off
            title(['Interface ', num2str(m)])
            xlabel('Position, x (nm)')
            ylabel('Carrier density (cm-3)')
            xlim([1e7*par.dcum(loc(m)-1)-0.5,1e7*par.dcum(loc(m))+0.5])
            ymin = min(min(min(n(:, pcum1(loc(m)-1):pcum1(loc(m)))), min(p(:, pcum1(loc(m)-1):pcum1(loc(m))))));
            ymax = max(max(max(n(:, pcum1(loc(m)-1):pcum1(loc(m)))), max(p(:, pcum1(loc(m)-1):pcum1(loc(m))))));
            ylim([ymin*0.1, ymax*10])
            
            %% Plot fluxes
            figure(601)
            subplot(1,length(loc),m)
            for i = 1:length(tarr)
                k = find(sol.t <= tarr(i));
                k = k(end);
                plot(x_sub*1e7, jn(k,:), x_sub*1e7, jp(k,:),...
                    x_sub*1e7, jn_ana(k,:), 'k-.', x_sub*1e7, jp_ana(k,:), 'k--')
                hold on
            end
            hold off
            title(['Interface ', num2str(m)])
            xlabel('Position, x (nm)')
            ylabel('Flux (cm-2s-1)')
            xlim([1e7*par.dcum(loc(m)-1)-0.5,1e7*par.dcum(loc(m))+0.5])
            
        end
    end
end
end