function compare_carrier_interfaces_6(sol, tarr)
% Function to compare the interfacial carrier densities between simulation
% and analytical solution
% v5 Calculated on the whole mesh, uses jn instead of js
% substitutes r for r=(jsp-jn)/x
% v6 Tests expressions for ns, ps for majority carrier sides of the
% interface, using xprime_n and xprime_p
% TARR = array of time points

par = sol.par;
dev = par.dev;
u = sol.u;
t = sol.t;
x = sol.x;
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
jn = zeros(1, length(x));
jp = zeros(1, length(x));
for i = 1:length(t)
    jn(i,:) = interp1(x_ihalf, j.n(i,:), x);
    jp(i,:) = interp1(x_ihalf, j.p(i,:), x);
end

% recombination
r_struct = dfana.calcr(sol);

ns1 = zeros(1, length(x));
ns2 = zeros(1, length(x));
ps1 = zeros(1, length(x));
ps2 = zeros(1, length(x));
jns1 = zeros(1, length(x));
jns2 = zeros(1, length(x));
jps1 = zeros(1, length(x));
jps2 = zeros(1, length(x));

for i = 1:length(tarr)
    pnt = find(sol.t <= tarr(i));
    pnt = pnt(end);
    
    for k = 1:length(x)
    alpha = alphaM(pnt,:);
    beta = betaM(pnt,:);
    
    %% Carrier densities
    if alpha(k) < 0 
        ns1(k) = u(pnt, par.pcum0(al-1)+1, 2);
        ns2(k) = u(pnt, par.pcum0(al+1)+1, 2);
    elseif alpha(k) >= 0
        ns1(k) = u(pnt, par.pcum0(al)+1, 2);
        ns2(k) = u(pnt, par.pcum0(al+2)+1, 2);
    end
    
    if beta(k) < 0      
        ps1(k) = u(pnt, par.pcum0(al-1)+1, 3);
        ps2(k) = u(pnt, par.pcum0(al+1)+1, 3);
    elseif beta(k) >= 0
        ps1(k) = u(pnt, par.pcum0(al)+1, 3);
        ps2(k) = u(pnt, par.pcum0(al+2)+1, 3);      
    end
    
    % Fluxes - approximate as calculated on half mesh
    if alpha(k) < 0
        jns1(k) = jn(pnt, par.pcum0(al-1)+1);
        jns2(k) = jn(pnt, par.pcum0(al+1)+1);
    elseif alpha(k) >= 0
        jns1(k) = jn(pnt, par.pcum0(al)+1);
        jns2(k) = jn(pnt, par.pcum0(al+2)+1);
    end
    
    if beta(k) < 0
        jps1(k) = jp(pnt, par.pcum0(al-1)+1);
        jps2(k) = jp(pnt, par.pcum0(al+1)+1);
    elseif beta(k) >= 0
        jps1(k) = jp(pnt, par.pcum0(al)+1);
        jps2(k) = jp(pnt, par.pcum0(al+2)+1);   
    end
    
    %recombination
    r = r_struct.tot(pnt, :);
    
    X = exp(-abs(alpha(k)).*xprime_n(k));
    Y = exp(-abs(beta(k)).*xprime_p(k));
    
    %% Carrier densities
    n1_ana(i,k) = ns1(k).*X - sign(alpha(k)).*(jn(pnt,k)./(-abs(alpha(k)).*mue(k).*Vt)).*(1 - X)...
        + sign(alpha(k)).*((jns1(k)-jn(pnt,k))./(alpha(k).^2.*mue(k).*Vt.*xprime_n(k))).*(1- X -abs(alpha(k)).*xprime_n(k).*X);
    
    n2_ana(i,k) = ns2(k).*X - sign(alpha(k)).*(jn(pnt,k)./(-abs(alpha(k)).*mue(k).*Vt)).*(1 - X)...
        + sign(alpha(k)).*((jns2(k)-jn(pnt,k))./(alpha(k).^2.*mue(k).*Vt.*xprime_n(k))).*(1- X -abs(alpha(k)).*xprime_n(k).*X);
    
    p1_ana(i,k) = ps1(k).*Y - sign(beta(k)).*(jp(pnt,k)./(-abs(beta(k)).*muh(k).*Vt)).*(1 - Y)...
        + sign(beta(k)).*((jps1(k)-jp(pnt,k))./(beta(k).^2.*muh(k).*Vt.*xprime_p(k))).*(1 - Y -abs(beta(k)).*xprime_p(k).*Y);
    
    p2_ana(i,k) = ps2(k).*Y - sign(beta(k)).*(jp(pnt,k)./(-abs(beta(k)).*muh(k).*Vt)).*(1 - Y)...
        + sign(beta(k)).*((jps2(k)-jp(pnt,k))./(beta(k).^2.*muh(k).*Vt.*xprime_p(k))).*(1 - Y -abs(beta(k)).*xprime_p(k).*Y);
    
    %% Components
    n1_comp1(i,k) = ns1(k).*X;
    n1_comp2(i,k) = - sign(alpha(k)).*(jn(pnt,k)./(-abs(alpha(k)).*mue(k).*Vt)).*(1 - X);
    n1_comp3(i,k) = + sign(alpha(k)).*((jns1(k)-jn(pnt,k))./(alpha(k).^2.*mue(k).*Vt.*xprime_n(k))).*(1- X -abs(alpha(k)).*xprime_n(k).*X);
    
    %% Fluxes
    jn1_ana(i,k) = jns1(k) + sign(alpha(k)).*r(k)*xprime_n(k);
    jn2_ana(i,k) = jns2(k) + sign(alpha(k)).*r(k)*xprime_n(k);
    jp1_ana(i,k) = jps1(k) + sign(beta(k)).*r(k)*xprime_p(k);
    jp2_ana(i,k) = jps2(k) + sign(beta(k)).*r(k)*xprime_p(k);
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
    semilogy(x*1e7, u(pnt, :, 2), x*1e7, u(pnt, :, 3),...
        x*1e7, n1_ana(i,:), 'k-.', x*1e7, p1_ana(i,:), 'k--')
    hold on
end
hold off
legend('n', 'p', 'n-ana', 'p-ana')
xlim([199.5,202.5])

% components
figure(5011)
for i = 1:length(tarr)
    pnt = find(sol.t <= tarr(i));
    pnt = pnt(end);
    plot(x*1e7, n1_comp1(i,:), x*1e7, n1_comp2(i,:), x*1e7, n1_comp3(i,:))
    hold on
end

%% Interface 2
figure(502)
for i = 1:length(tarr)
    pnt = find(sol.t <= tarr(i));
    pnt = pnt(end);
    semilogy(x*1e7, u(pnt, :, 2), x*1e7, u(pnt, :, 3),...
        x*1e7, n2_ana(i,:), 'k-.', x*1e7, p2_ana(i,:), 'k--')
    hold on
end
hold off
legend('n', 'p', 'n-ana', 'p-ana')
xlim([541.5,544.5])

%% Plot fluxes
%% Interface 1
figure(503)

for i = 1:length(tarr)
    pnt = find(sol.t <= tarr(i));
    pnt = pnt(end);

    plot(x*1e7, jn(pnt,:), x*1e7, jp(pnt,:),... 
        x*1e7, jn1_ana(i,:), 'k-.', x*1e7, jp1_ana(i,:), 'k--')
    hold on
end
xlim([199.5,202.5])
hold off
   
%% Interface 2
figure(504)
hold on
for i = 1:length(tarr)
    pnt = find(sol.t <= tarr(i));
    pnt = pnt(end);

    plot(x*1e7, jn(pnt,:), x*1e7, jp(pnt,:),... 
         x*1e7, jn2_ana(i,:), 'k-.', x*1e7, jp2_ana(i,:), 'k--')
end
xlim([541.5,544.5])
hold off
end



