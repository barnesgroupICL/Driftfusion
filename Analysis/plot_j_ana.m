function plot_j_ana(sol_df, tpoint)
% Only works for HTL-i-ETL structure

% Script to compare the interfacial recombination fluxes from DF and IM
par = sol_df.par;
u = sol_df.u;
t = sol_df.t;
x = sol_df.x;
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
alpha = par.dev.alpha;
beta = par.dev.beta;
xprime_n = par.dev.xprime_n;
xprime_p = par.dev.xprime_p;
mue = par.dev.mue;
muh = par.dev.muh;
kBT = par.kB*par.T;
x = sol_df.x;

%% Driftfusion 2D
ns1_df = u(:, par.pcum(al-1)+1, 2);
ps1_df = u(:, par.pcum(al-2)+1, 3);
ns2_df = u(:, par.pcum(al+1)+1, 2);
ps2_df = u(:, par.pcum(al)+1, 3);

n_ana = ns1_df.*exp(-abs(alpha).*xprime_n);
jn_ana_drift = mue.*kBT.*abs(alpha).*n_ana;
jn_ana_diff = -mue.*kBT.*abs(alpha).*n_ana;

for i=1:length(sol_df.t)
    dn_anadx(i,:) = gradient(n_ana(i,:), xprime_n);
    jn_ana(i,:) = mue.*kBT.*(abs(alpha).*n_ana(i,:) - dn_anadx(i,:));
end

p_ana = ps1_df.*exp(-abs(beta).*xprime_p);
jp_ana = 2*mue.*kBT.*abs(beta).*p_ana;

% find the time
p1 = find(sol_df.t <= tpoint);
p1 = p1(end);

%% Plot the carrier densities
dfplot.npx(sol_df, tpoint)
hold on
plot(x*1e7, n_ana(p1,:), '-.', x*1e7, p_ana(p1,:), '-.')
xlim([x(par.pcum(al-2))*1e7-1, x(par.pcum(al-1))*1e7+1])
hold off

%% Plot the fluxes
dfplot.jx(sol_df, tpoint)
hold on
plot(x*1e7, jn_ana(p1,:),'--');%, x*1e7, jp_ana(p1,:))
%xlim([x(par.pcum(al-2))*1e7-1, x(par.pcum(al-1))*1e7+1])
hold off

%%Plot DD currents
dfplot.Jddx(sol_df, tpoint)
hold on
plot(x*1e7, -par.e*jn_ana_drift(p1,:),'--', x*1e7, -par.e*jn_ana_diff(p1,:));%, x*1e7, jp_ana(p1,:))
%xlim([x(par.pcum(al-2))*1e7-1, x(par.pcum(al-1))*1e7+1])
hold off

%%Plot DD currents
dfplot.Jx(sol_df, tpoint)
end

