function compare_carrier_interfaces_9(sol, tarr)
% Function to compare the interfacial carrier densities between simulation
% and analytical solution
% v5 Calculated on the whole mesh, uses jn instead of js
% substitutes r for r=(jsp-jn)/x_ihalf
% v6 Tests expressions for ns, ps for majority carrier sides of the
% v7 as with v6 but calculated on the half mesh
% v7.1 as with v7 but recombination flux terms are taken from LHS only
% interface, using xprime_n and xprime_p
% v9 As with v7.1 but attempt to remove the sign functions and simply use
% conditional statement to switch between expressions.
% TARR = array of time points

par = sol.par;
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
dint = dev.dint;
mue = dev.mue;
muh = dev.muh;
Vt = par.kB*par.T;
pcum0 = par.pcum0;

for i=1:length(t)
    [~, dVdx(i,:)] = pdeval(0, x, sol.u(i,:,1), x_ihalf);
    [n(i,:), ~] = pdeval(0, x, sol.u(i,:,2), x_ihalf);
    [p(i,:), ~] = pdeval(0, x, sol.u(i,:,3), x_ihalf);
end

alphaM = par.q*dVdx/(par.kB*par.T) + alpha_prime;
betaM = par.q*-dVdx/(par.kB*par.T) + beta_prime; 

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

ns1 = zeros(1, length(x_ihalf));
ns2 = zeros(1, length(x_ihalf));
ps1 = zeros(1, length(x_ihalf));
ps2 = zeros(1, length(x_ihalf));
jns1 = zeros(1, length(x_ihalf));
jns2 = zeros(1, length(x_ihalf));
jps1 = zeros(1, length(x_ihalf));
jps2 = zeros(1, length(x_ihalf));

for i = 1:length(tarr)
    pnt = find(sol.t <= tarr(i));
    pnt = pnt(end);
    
    for k = 1:length(x_ihalf)
    alpha = alphaM(pnt,:);
    beta = betaM(pnt,:);
    
    jns1 = jn(pnt, par.pcum0(al-1)+1);
    jns2 = jn(pnt, par.pcum0(al+1)+1);
    jps1 = jp(pnt, par.pcum0(al-1)+1);
    jps2 = jp(pnt, par.pcum0(al+1)+1);
        
    %recombination
    r = r_struct.tot(pnt, :);
    %% Carrier densities
    if alpha(k) < 0 
        ns1 = u(pnt, par.pcum0(al-1)+1, 2);
        ns2 = u(pnt, par.pcum0(al+1)+1, 2);
    
        X = exp(alpha(k)*xprime(k));    
        n1_ana(i,k) = ns1*X + (jn(pnt,k)/(alpha(k)*mue(k)*Vt))*(1 - X)...
            - ((jns1-jn(pnt,k))/(alpha(k)^2*mue(k)*Vt*xprime(k)))*(1- X + alpha(k)*xprime(k)*X);
    
        n2_ana(i,k) = ns2*X + (jn(pnt,k)/(alpha(k)*mue(k)*Vt))*(1 - X)...
            - ((jns2-jn(pnt,k))/(alpha(k)^2*mue(k)*Vt*xprime(k)))*(1- X + alpha(k)*xprime(k)*X);

        %% Fluxes
        jn1_ana(i,k) = jns1 - r(k)*xprime_n(k);
        jn2_ana(i,k) = jns2 - r(k)*xprime_n(k);
        
    elseif alpha(k) >= 0
        X = exp(alpha(k)*xprime(k));  
        ns1 = u(pnt, par.pcum0(al)+1, 2);
        ns2 = u(pnt, par.pcum0(al+2)+1, 2);
        
        grad_jp1 = (jps1-jp(pnt,k))/xprime_p(k);
        grad_jp2 = (jps2-jp(pnt,k))/xprime_p(k);
        
        n1_ana(i,k) = NaN;
        n2_ana(i,k) = NaN;
        
        %% Fluxes
        jn1_ana(i,k) = jns1 + r(k)*xprime_n(k);
        jn2_ana(i,k) = jns2 + r(k)*xprime_n(k);
    else
        X = exp(alpha(k)*xprime(k));  
        n1_ana(i,k) = NaN;
        n2_ana(i,k) = NaN;   
    end
    
    if beta(k) < 0      
        ps1 = u(pnt, par.pcum0(al-1)+1, 3);
        ps2 = u(pnt, par.pcum0(al+1)+1, 3);
        Y = exp(beta(k)*xprime_p(k));
        
        p1_ana(i,k) = NaN;
        p2_ana(i,k) = NaN;
        
        jp1_ana(i,k) = jps1 - r(k)*xprime_p(k);
        jp2_ana(i,k) = jps2 - r(k)*xprime_p(k);
    elseif beta(k) >= 0
        ps1 = u(pnt, par.pcum0(al)+1, 3);
        ps2 = u(pnt, par.pcum0(al+2)+1, 3);
        
        grad_jn1 = (jns1-jn(pnt,k))/xprime_n(k);
        grad_jn2 = (jns2-jn(pnt,k))/xprime_n(k);
        
        grad_jp1 = (jps1-jp(pnt,k))/xprime_p(k);
        grad_jp2 = (jps2-jp(pnt,k))/xprime_p(k);
        
        Y = exp(-beta(k)*(dint(k)-xprime(k)));
        
        p1_ana(i,k) = ps1*Y + (jp(pnt,k)/(beta(k).*muh(k).*Vt))*(1 - Y)...
        	+ (grad_jp1/(beta(k)^2*muh(k)*Vt))*(1 - Y - beta(k)*(dint(k)-xprime(k))*Y);
    
        p2_ana(i,k) = ps2*Y + (jp(pnt,k)/(beta(k)*muh(k)*Vt))*(1 - Y)...
            + (grad_jn2/(beta(k)^2*muh(k)*Vt))*(1 - Y - beta(k)*(dint(k)-xprime(k))*Y);
        
        jp1_ana(i,k) = jps1 + r(k)*xprime_p(k);
        jp2_ana(i,k) = jps2 + r(k)*xprime_p(k);
    else
        Y = exp(beta(k)*xprime_p(k));
        
        p1_ana(i,k) = NaN;
        p2_ana(i,k) = NaN;
    end

    %% Components
    n1_comp1(i,k) = ns1*X;
    n1_comp2(i,k) = - sign(alpha(k))*(jn(pnt,k)/(-abs(alpha(k))*mue(k)*Vt))*(1 - X);
    n1_comp3(i,k) = + sign(alpha(k))*((jns1-jn(pnt,k))/(alpha(k).^2*mue(k)*Vt*xprime_n(k)))*(1- X -abs(alpha(k))*xprime_n(k)*X);
 
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

%% Plot carrier densities
%% Interface 1
figure(501)
for i = 1:length(tarr)
    pnt = find(sol.t <= tarr(i));
    pnt = pnt(end);
    subplot(2,2,1)
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
    subplot(2,2,2)
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

%% Plot fluxes
%% Interface 1
for i = 1:length(tarr)
    pnt = find(sol.t <= tarr(i));
    pnt = pnt(end);
    subplot(2,2,3)
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
    subplot(2,2,4)
    plot(x_ihalf*1e7, jn(pnt,:), x_ihalf*1e7, jp(pnt,:),... 
         x_ihalf*1e7, jn2_ana(i,:), 'k-.', x_ihalf*1e7, jp2_ana(i,:), 'k--')
    hold on
end
title('Interface 2')
xlabel('Position, x (nm)')
ylabel('Flux (cm-2s-1)')
xlim([1e7*par.dcum0(al+1)-0.5,1e7*par.dcum0(al+2)+0.5])
hold off

% components
figure(5011)
for i = 1:length(tarr)
    pnt = find(sol.t <= tarr(i));
    pnt = pnt(end);
    plot(x_ihalf*1e7, n1_comp1(i,:), x_ihalf*1e7, n1_comp2(i,:), x_ihalf*1e7, n1_comp3(i,:), x_ihalf*1e7, n1_ana(i,:), 'k-.')
    hold on
end
xlabel('Position, x (nm)')
ylabel('Carrier density (cm-3)')
legend('n comp', 'j comp', 'r comp', 'sum')
xlim([1e7*par.dcum0(al-1)-0.5,1e7*par.dcum0(al)+0.5])
hold off
end



