function compare_carrier_interfaces_5(sol, tarr)
% Function to compare the interfacial carrier densities between simulation
% and analytical solution
% v5 Calculated on the whole mesh, uses jn instead of js
% substitutes r for r=(jsp-jn)/x
% TARR = array of time points

par = sol.par;
dev = par.dev;
u = sol.u;
t = sol.t;
x = sol.x;
al = par.active_layer;
ni = par.ni(al);
sn1 = par.sn(al-1);
sp1 = par.sp(al-1);
sn2 = par.sn(al+1);
sp2 = par.sp(al+1);
nt1 = par.nt(al-1);
pt1 = par.pt(al-1);
nt2 = par.nt(al+1);
pt2 = par.pt(al+1);
xprime_p = dev.xprime;
xprime_n = dev.xprime;
alpha_prime = dev.alpha_prime;
beta_prime = dev.beta_prime;
mue = dev.mue;
muh = dev.muh;
Vt = par.kB*par.T;
pcum0 = par.pcum0;
for i=1:length(t)
    [~, dVdx(i,:)] = pdeval(0, x, sol.u(i,:,1), x);
end

alphaM = par.q*dVdx./(par.kB*par.T) + alpha_prime;
betaM = par.q*-dVdx./(par.kB*par.T) + beta_prime; 

% Fluxes - approximate as on half mesh
[J, j, x_ihalf] = dfana.calcJ(sol);
for i = 1:length(t)
    jn(i,:) = interp1(x_ihalf, j.n(i,:), x);
    jp(i,:) = interp1(x_ihalf, j.p(i,:), x);
end

% recombination
r_struct = dfana.calcr(sol);
    
for i = 1:length(tarr)
    p1 = find(sol.t <= tarr(i));
    p1 = p1(end);

    alpha = alphaM(p1,:);
    beta = betaM(p1,:);
    
    %% Carrier densities
    ns1 = u(p1, par.pcum0(al-1)+1, 2);
    ps1 = u(p1, par.pcum0(al-1)+1, 3);
    ns2 = u(p1, par.pcum0(al+1)+1, 2);
    ps2 = u(p1, par.pcum0(al+1)+1, 3);

    % Fluxes - approximate as calculated on half mesh
    jns1 = jn(p1, par.pcum0(al-1)+2);
    jps1 = jp(p1, par.pcum0(al-1)+2);
    jns2 = jn(p1, par.pcum0(al+1)+2);
    jps2 = jp(p1, par.pcum0(al+1)+2);  
    
    %recombination
    r = r_struct.tot(p1, :);
    %r = 0;
    
    ps1_ana(i,:) = ps1.*exp(beta.*xprime_p) + jp(p1,:)./(beta.*muh.*Vt).*(1-exp(beta.*xprime_p))...
        - ((jps1-jp(p1,:))./(beta.^2.*muh.*Vt.*xprime_p)).*(1-exp(beta.*xprime_p)+ beta.*xprime_p.*exp(beta.*xprime_p));
    
    ps2_ana(i,:) = ps2.*exp(beta.*xprime_p) + jp(p1,:)./(beta.*muh.*Vt).*(1-exp(beta.*xprime_p))...
        - ((jps2-jp(p1,:))./(beta.^2.*muh.*Vt.*xprime_p)).*(1-exp(beta.*xprime_p)+ beta.*xprime_p.*exp(beta.*xprime_p));
    
    ns1_ana(i,:) = ns1.*exp(alpha.*xprime_n) + jn(p1,:)./(alpha.*mue.*Vt).*(1-exp(alpha.*xprime_n))...
        - ((jns1-jn(p1,:))./(alpha.^2.*mue.*Vt.*xprime_n)).*(1-exp(alpha.*xprime_n)+ alpha.*xprime_n.*exp(alpha.*xprime_n));
   
    
    ns2_ana(i,:) = ns2.*exp(alpha.*xprime_n) + jn(p1,:)./(alpha.*mue.*Vt).*(1-exp(alpha.*xprime_n))...
        - ((jns2-jn(p1,:))./(alpha.^2.*mue.*Vt.*xprime_n)).*(1-exp(alpha.*xprime_n)+ alpha.*xprime_n.*exp(alpha.*xprime_n));
end

%clean up for plotting
ps1_ana(:, pcum0(1):pcum0(2)-1) = NaN;
ns1_ana(:, pcum0(1):pcum0(2)-1) = NaN;
ps1_ana(:, pcum0(3)+1:end) = NaN;
ns1_ana(:, pcum0(3)+1:end) = NaN;
ps2_ana(:, pcum0(1):pcum0(4)-1) = NaN;
ns2_ana(:, pcum0(1):pcum0(4)-1) = NaN;
ps2_ana(:, pcum0(5)+1:end) = NaN;
ns2_ana(:, pcum0(5)+1:end) = NaN;

% plot carrier densities
dfplot.npx(sol, tarr)
hold on
for i = 1:length(tarr)
    plot(x*1e7, ps1_ana(i,:), 'k--', x*1e7, ps2_ana(i,:), 'k--',...
        x*1e7, ns1_ana(i,:), 'k-.', x*1e7, ns2_ana(i,:), 'k-.')
end
% legend('HTL', 'Int 1', 'AL', 'Int 2', 'ETL',...
%     'n', 'p', 'ps1-ana', 'ps2-ana', 'ns1-ana', 'ns2-ana')
xlim([199.5,202.5])
xlim([541.5,544.5])

hold off
    
end



