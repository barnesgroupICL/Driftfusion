function compare_carrier_interfaces_8_2(sol, tarr)
% Function to compare the interfacial carrier densities between simulation
% and analytical solution
% v5 Calculated on the whole mesh, uses jn instead of js
% substitutes r for r=(jsp-jn)/x_ihalf
% v6 Tests expressions for ns, ps for majority carrier sides of the
% v7 as with v6 but calculated on the half mesh
% v8 uses left hand boundary values only and adapts expression for n(d) =
% ns and p(d) = ps dependent on sign of alpha and beta
% v8.2 includes a switch to plot ns ps calculated from n(x), p(x), jn(x),
% jp(x), jns and jps
% interface, using xprime_n and xprime_p
% TARR = array of time points
plotswitch = 2;     % 1 = plot n(x), p(x), 2 = plot ns, ps

par = sol.par;
T = par.T;
kB = par.kB;
dev = par.dev_ihalf;
u = sol.u;
t = sol.t;
x = sol.x;
x_ihalf = par.x_ihalf;
al = par.active_layer;
ni = par.ni(al);
sn1 = par.sn(al-1);
spnt = par.sp(al-1);
sn2 = par.sn(al+1);
sp2 = par.sp(al+1);
nt1 = par.nt(al-1);
pt1 = par.pt(al-1);
nt2 = par.nt(al+1);
pt2 = par.pt(al+1);
xprime_p = dev.xprime_p;
xprime_n = dev.xprime_n;
xprime = dev.xprime;
alpha_prime = dev.alpha_prime;
beta_prime = dev.beta_prime;
mue = dev.mue;
muh = dev.muh;
dint = dev.dint;
Vt = par.kB*par.T;
kB = par.kB;
T = par.T;
pcum0 = par.pcum0;

for i=1:length(t)
    [~, dVdx(i,:)] = pdeval(0, x, sol.u(i,:,1), x_ihalf);
    [n(i,:), ~] = pdeval(0, x, sol.u(i,:,2), x_ihalf);
    [p(i,:), ~] = pdeval(0, x, sol.u(i,:,3), x_ihalf);
end

alphaM = par.q*dVdx./(par.kB*par.T) + alpha_prime;
betaM = par.q*-dVdx./(par.kB*par.T) + beta_prime;

% Fluxes - approximate as on half mesh
[J, j, x_ihalf] = dfana.calcJ(sol);
jn = zeros(1, length(x_ihalf));
jp = zeros(1, length(x_ihalf));
for i = 1:length(t)
    jn(i,:) = j.n(i,:);
    jp(i,:) = j.p(i,:);
end

% recombination
r_struct = dfana.calcr(sol);

for i = 1:length(tarr)
    pnt = find(sol.t <= tarr(i));
    pnt = pnt(end);
    
    for k = 1:length(x_ihalf)
        alpha = alphaM(pnt,:);
        beta = betaM(pnt,:);
        
        %% Carrier densities
        if alpha(k) < 0
            ns1 = u(pnt, par.pcum0(al-1), 2);
            ns2 = u(pnt, par.pcum0(al+1), 2);
        elseif alpha(k) >= 0
            ns1 = u(pnt, par.pcum0(al)+1, 2);
            ns2 = u(pnt, par.pcum0(al+2)+1, 2);
        end
        
        if beta(k) < 0
            ps1 = u(pnt, par.pcum0(al-1), 3);
            ps2 = u(pnt, par.pcum0(al+1), 3);
        elseif beta(k) >= 0
            ps1 = u(pnt, par.pcum0(al)+1, 3);
            ps2 = u(pnt, par.pcum0(al+2)+1, 3);
        end
        
        % Fluxes - approximate as calculated on half mesh
        %     if alpha(k) < 0
        jns1 = jn(pnt, par.pcum0(al-1)+1);
        jns2 = jn(pnt, par.pcum0(al+1)+1);
        %     elseif alpha(k) >= 0
        %         jns1 = jn(pnt, par.pcum0(al)-1);
        %         jns2 = jn(pnt, par.pcum0(al+2)-1);
        %     end
        
        %     if beta(k) < 0
        jps1 = jp(pnt, par.pcum0(al-1)+1);
        jps2 = jp(pnt, par.pcum0(al+1)+1);
        %     elseif beta(k) >= 0
        %         jps1 = jp(pnt, par.pcum0(al)-1);
        %         jps2 = jp(pnt, par.pcum0(al+2)-1);
        %     end
        
        %recombination
        r = r_struct.tot(pnt, :);
        
        if plotswitch == 1
            
            %% Carrier densities
            if alpha(k) <0
                X = exp(alpha(k).*xprime(k));
                n1_ana(i,k) = ns1*X + jns1/(alpha(k)*mue(k)*Vt)*(1 - X)...
                    - ((jns1-jn(pnt,k))/(alpha(k)^2*mue(k)*Vt*xprime(k)))*(1 - X + alpha(k)*xprime(k));
                
                n2_ana(i,k) = ns2*X + jns2/(alpha(k)*mue(k)*Vt)*(1 - X)...
                    - ((jns2-jn(pnt,k))/(alpha(k)^2*mue(k)*Vt*xprime(k)))*(1 - X + alpha(k)*xprime(k));
                
            elseif alpha(k) > 0
                Xstar = exp(alpha(k)*(xprime(k)-dint(k)));
                n1_ana(i,k) = ns1*Xstar + (jns1/(alpha(k)*mue(k)*Vt))*(1 - Xstar)...
                    - ((jns1-jn(pnt,k))/(alpha(k)^2*mue(k)*Vt*xprime(k)))*(1 + alpha(k)*xprime(k) - (1 + alpha(k)*dint(k))*Xstar);
                
                n2_ana(i,k) = ns2*Xstar + (jns2/(alpha(k)*mue(k)*Vt))*(1 - Xstar)...
                    - ((jns2-jn(pnt,k))/(alpha(k)^2*mue(k)*Vt*xprime(k)))*(1 + alpha(k)*xprime(k) - (1 + alpha(k)*dint(k))*Xstar);
            end
            
            if beta(k) <0
                Y = exp(beta(k)*xprime(k));
                p1_ana(i,k) = ps1*Y + (jps1/(beta(k)*muh(k)*Vt))*(1 - Y)...
                    - ((jps1-jp(pnt,k))/(beta(k)^2*muh(k)*Vt*xprime(k)))*(1 - Y + beta(k)*xprime(k));
                
                p2_ana(i,k) = ps2*Y + (jps2/(beta(k)*muh(k)*Vt))*(1 - Y)...
                    - ((jps2-jp(pnt,k))/(beta(k)^2*muh(k)*Vt*xprime(k)))*(1 - Y + beta(k)*xprime(k));
                
            elseif beta(k) >= 0
                Ystar = exp(beta(k)*(xprime(k)-dint(k)));
                p1_ana(i,k) = ps1*Ystar + (jps1/(beta(k)*muh(k)*Vt))*(1 - Ystar)...
                    - ((jps1-jp(pnt,k))/(beta(k)^2*muh(k)*Vt*xprime(k)))*(1 + beta(k)*xprime(k) - (1 + beta(k)*dint(k))*Ystar);
                
                p2_ana(i,k) = ps2*Ystar + (jps2/(beta(k)*muh(k)*Vt))*(1 - Ystar)...
                    - ((jps2-jp(pnt,k))/(beta(k)^2*muh(k)*Vt*xprime(k)))*(1 + beta(k)*xprime(k) - (1 + beta(k)*dint(k))*Ystar);
            end
            
        end
        
        if plotswitch == 2
            %% Carrier densities
%             if alpha(k) < 0
                X = exp(alpha(k).*xprime(k));
                ns1_ana_l(i,k) = (1/X)*(n(pnt,k) - (jns1/(alpha(k)*mue(k)*Vt))*(1 - X) + ((jns1-jn(pnt,k))/(alpha(k)^2*mue(k)*Vt*xprime(k)))*(1 - X + alpha(k)*xprime(k)));
                ns2_ana_l(i,k) = (1/X)*(n(pnt,k) - (jns2/(alpha(k)*mue(k)*Vt))*(1 - X) + ((jns2-jn(pnt,k))/(alpha(k)^2*mue(k)*Vt*xprime(k)))*(1 - X + alpha(k)*xprime(k)));
                
%             elseif alpha(k) >= 0
                Xstar = exp(alpha(k)*(xprime(k)-dint(k)));
                ns1_ana_r(i,k) = (1/Xstar)*(n(pnt,k) - (jns1/(alpha(k)*mue(k)*Vt))*(1 - Xstar) + ((jns1-jn(pnt,k))/(alpha(k)^2*mue(k)*Vt*xprime(k)))*(1 + alpha(k)*xprime(k) - (1 + alpha(k)*dint(k))*Xstar));
                ns2_ana_r(i,k) = (1/Xstar)*(n(pnt,k) - (jns2/(alpha(k)*mue(k)*Vt))*(1 - Xstar) + ((jns2-jn(pnt,k))/(alpha(k)^2*mue(k)*Vt*xprime(k)))*(1 + alpha(k)*xprime(k) - (1 + alpha(k)*dint(k))*Xstar));
%             end
                ns1_ana(i,k) = max(ns1_ana_l(i,k), ns1_ana_r(i,k));
                ns2_ana(i,k) = max(ns2_ana_l(i,k), ns2_ana_r(i,k));
                
%             if beta(k) <0
                Y = exp(beta(k)*xprime(k));
                ps1_ana_l(i,k) = (1/Y)*(p(pnt,k) - (jps1/(beta(k)*muh(k)*Vt))*(1 - Y) + ((jps1-jp(pnt,k))/(beta(k)^2*muh(k)*Vt*xprime(k)))*(1 - Y + beta(k)*xprime(k)));
                ps2_ana_l(i,k) = (1/Y)*(p(pnt,k) - (jps2/(beta(k)*muh(k)*Vt))*(1 - Y) + ((jps2-jp(pnt,k))/(beta(k)^2*muh(k)*Vt*xprime(k)))*(1 - Y + beta(k)*xprime(k)));

%             elseif beta(k) >= 0
                Ystar = exp(beta(k)*(xprime(k)-dint(k)));
                ps1_ana_r(i,k) = (1/Ystar)*(p(pnt,k) - (jps1/(beta(k)*muh(k)*Vt))*(1 - Ystar) + ((jps1-jp(pnt,k))/(beta(k)^2*muh(k)*Vt*xprime(k)))*(1 + beta(k)*xprime(k) - (1 + beta(k)*dint(k))*Ystar));
                ps2_ana_r(i,k) = (1/Ystar)*(p(pnt,k) - (jps2/(beta(k)*muh(k)*Vt))*(1 - Ystar) + ((jps2-jp(pnt,k))/(beta(k)^2*muh(k)*Vt*xprime(k)))*(1 + beta(k)*xprime(k) - (1 + beta(k)*dint(k))*Ystar));
%             end
                ps1_ana(i,k) = max(ps1_ana_l(i,k), ps1_ana_r(i,k));
                ps2_ana(i,k) = max(ps2_ana_l(i,k), ps2_ana_r(i,k));   
        end
        
        %% Fluxes
        jn1_ana(i,k) = jns1 - r(k)*xprime(k);
        jn2_ana(i,k) = jns2 - r(k)*xprime(k);
        jp1_ana(i,k) = jps1 - r(k)*xprime(k);
        jp2_ana(i,k) = jps2 - r(k)*xprime(k);
    end
end

%clean up for plotting
p1_ana(:, pcum0(1):pcum0(2)) = NaN;
n1_ana(:, pcum0(1):pcum0(2)) = NaN;
p1_ana(:, pcum0(3)+1:end) = NaN;
n1_ana(:, pcum0(3)+1:end) = NaN;
p2_ana(:, pcum0(1):pcum0(4)) = NaN;
n2_ana(:, pcum0(1):pcum0(4)) = NaN;
p2_ana(:, pcum0(5)+1:end) = NaN;
n2_ana(:, pcum0(5)+1:end) = NaN;

if plotswitch == 1
    %% Plot carrier densities
    %% Interface 1
    figure(501)
    for i = 1:length(tarr)
        pnt = find(sol.t <= tarr(i));
        pnt = pnt(end);
        subplot(1,2,1)
        semilogy(x_ihalf*1e7, n(pnt, :), x_ihalf*1e7, p(pnt, :),...
            x_ihalf*1e7, n1_ana(i,:), 'k-.', x_ihalf*1e7, p1_ana(i,:), 'k--')
        hold on
    end
    hold off
    title('Interface 1')
    %legend('n', 'p', 'n-ana', 'p-ana')
    xlabel('Position, x (nm)')
    ylabel('Carrier density (cm-3)')
    xlim([1e7*par.dcum0(al-1)-0.5,1e7*par.dcum0(al)+0.5])
    
    %% Interface 2
    for i = 1:length(tarr)
        pnt = find(sol.t <= tarr(i));
        pnt = pnt(end);
        subplot(1,2,2)
        semilogy(x_ihalf*1e7, n(pnt, :), x_ihalf*1e7, p(pnt, :),...
            x_ihalf*1e7, n2_ana(i,:), 'k-.', x_ihalf*1e7, p2_ana(i,:), 'k--')
        hold on
    end
    hold off
    title('Interface 2')
    xlabel('Position, x (nm)')
    ylabel('Carrier density (cm-3)')
    %legend('n', 'p', 'n-ana', 'p-ana')
    xlim([1e7*par.dcum0(al+1)-0.5,1e7*par.dcum0(al+2)+0.5])
    
    figure(502)
    %% Plot fluxes
    %% Interface 1
    for i = 1:length(tarr)
        pnt = find(sol.t <= tarr(i));
        pnt = pnt(end);
        subplot(1,2,1)
        plot(x_ihalf*1e7, jn(pnt,:), x_ihalf*1e7, jp(pnt,:),...
            x_ihalf*1e7, jn1_ana(i,:), 'k-.', x_ihalf*1e7, jp1_ana(i,:), 'k--')
        hold on
    end
    title('Interface 1')
    xlabel('Position, x (nm)')
    ylabel('Flux (cm-2s-1)')
    xlim([1e7*par.dcum0(al-1)-0.5,1e7*par.dcum0(al)+0.5])
    hold off
    
    %% Interface 2
    for i = 1:length(tarr)
        pnt = find(sol.t <= tarr(i));
        pnt = pnt(end);
        subplot(1,2,2)
        plot(x_ihalf*1e7, jn(pnt,:), x_ihalf*1e7, jp(pnt,:),...
            x_ihalf*1e7, jn2_ana(i,:), 'k-.', x_ihalf*1e7, jp2_ana(i,:), 'k--')
        hold on
    end
    title('Interface 2')
    xlabel('Position, x (nm)')
    ylabel('Flux (cm-2s-1)')
    xlim([1e7*par.dcum0(al+1)-0.5,1e7*par.dcum0(al+2)+0.5])
    hold off
    
end

if plotswitch ==2
    figure(502)
    for i = 1:length(tarr)
        pnt = find(sol.t <= tarr(i));
        pnt = pnt(end);
        subplot(1,2,1)
        semilogy(x_ihalf*1e7, n(pnt, :), x_ihalf*1e7, p(pnt, :),...
            x_ihalf*1e7, ns1_ana(i,:), 'k-.', x_ihalf*1e7, ps1_ana(i,:), 'k--')
        hold on
    end
    hold off
    title('Interface 1')
    %legend('n', 'p', 'n-ana', 'p-ana')
    xlabel('Position, x (nm)')
    ylabel('Carrier density (cm-3)')
    xlim([1e7*par.dcum0(al-1)-0.5,1e7*par.dcum0(al)+0.5])
    
    for i = 1:length(tarr)
        pnt = find(sol.t <= tarr(i));
        pnt = pnt(end);
        subplot(1,2,2)
        semilogy(x_ihalf*1e7, n(pnt, :), x_ihalf*1e7, p(pnt, :),...
            x_ihalf*1e7, ns2_ana(i,:), 'k-.', x_ihalf*1e7, ps2_ana(i,:), 'k--')
        hold on
    end
    hold off
    title('Interface 1')
    %legend('n', 'p', 'n-ana', 'p-ana')
    xlabel('Position, x (nm)')
    ylabel('Carrier density (cm-3)')
    xlim([1e7*par.dcum0(al+1)-0.5,1e7*par.dcum0(al+2)+0.5])
    
    % % components
    % figure(5011)
    % for i = 1:length(tarr)
    %     pnt = find(sol.t <= tarr(i));
    %     pnt = pnt(end);
    %     plot(x_ihalf*1e7, n1_comp1(i,:), x_ihalf*1e7, n1_comp2(i,:), x_ihalf*1e7, n1_comp3(i,:), x_ihalf*1e7, n1_ana(i,:), 'k-.')
    %     hold on
    % end
    % xlabel('Position, x (nm)')
    % ylabel('Carrier density (cm-3)')
    % legend('n comp', 'j comp', 'r comp', 'sum')
    % xlim([199.5,202.5])
    % hold off
end



