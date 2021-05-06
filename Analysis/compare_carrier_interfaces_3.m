function compare_carrier_interfaces_3(sol, tarr)
% Function to compare the interfacial carrier densities between simulation
% and analytical solution
% v3 calculates on the half mesh
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

for i=1:length(t)
    [n(i,:), ~] = pdeval(0, x, u(i,:,2), x_ihalf);
    [p(i,:), ~] = pdeval(0, x, u(i,:,3), x_ihalf);
    [~, dVdx(i,:)] = pdeval(0, x, sol.u(i,:,1), x_ihalf);
end

alpha = par.q*dVdx./(par.kB*par.T) + alpha_prime;
beta = par.q*-dVdx./(par.kB*par.T) + beta_prime; 

% Fluxes - approximate as on half mesh
[J, j, x_ihalf] = dfana.calcJ(sol);
% recombination
r_struct = dfana.calcr(sol);
    
for i = 1:length(tarr)
    p1 = find(sol.t <= tarr(i));
    p1 = p1(end);

    alpha = alpha(p1,:);
    beta = beta(p1,:);
    
    %% Carrier densities
    ns1 = n(p1, par.pcum0(al-1)+1);
    ps1 = p(p1, par.pcum0(al-1)+1);
    ns2 = n(p1, par.pcum0(al+1)+1);
    ps2 = p(p1, par.pcum0(al+1)+1);

    % Fluxes
    jns1 = j.n(p1, par.pcum0(al-1)+1);
    jps1 = j.p(p1, par.pcum0(al-1)+1);
    jns2 = j.n(p1, par.pcum0(al+1)+1);
    jps2 = j.p(p1, par.pcum0(al+1)+1);  
    
    %recombination
    r = r_struct.tot(p1, :);
    
    ps1_ana = ps1.*exp(beta.*xprime_p) + jps1./(beta.*muh.*Vt).*(1-exp(beta.*xprime_p))...
        - (r./(beta.^2.*muh.*Vt)).*(1-exp(beta.*xprime_p)+ beta.*xprime_p);
    
    ps2_ana = ps2.*exp(beta.*xprime_p) + jps2./(beta.*muh.*Vt).*(1-exp(beta.*xprime_p))...
        - (r./(beta.^2.*muh.*Vt)).*(1-exp(beta.*xprime_p)+ beta.*xprime_p);

    ns1_ana = ns1.*exp(alpha.*xprime_n) + jns1./(alpha.*mue.*Vt).*(1-exp(alpha.*xprime_n));...
        - (r./(alpha.^2.*mue.*Vt)).*(1-exp(alpha.*xprime_n)+ alpha.*xprime_n);
    
    ns2_ana = ns2.*exp(alpha.*xprime_n) + jns2./(alpha.*mue.*Vt).*(1-exp(alpha.*xprime_n))...
        - (r./(alpha.^2.*mue.*Vt)).*(1-exp(alpha.*xprime_n)+ alpha.*xprime_n);
end

%clean up for plotting
ps1_ana(par.pcum0(al)+1:end) = 0;
ns1_ana(par.pcum0(al)+1:end) = 0;
ps2_ana(par.pcum0(1):par.pcum0(al)+1) = 0;
ns2_ana(par.pcum0(1):par.pcum0(al)+1) = 0;
% plot carrier densities
dfplot.npx(sol, tarr)
hold on
for i = 1:length(tarr)
    plot(x_ihalf*1e7, ps1_ana, 'k--', x_ihalf*1e7, ps2_ana, 'k--',...
        x_ihalf*1e7, ns1_ana, 'k-.', x_ihalf*1e7, ns2_ana, 'k-.')
end
hold off
    
end



