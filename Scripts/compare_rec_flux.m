function compare_rec_flux(sol_df, sol_im)
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

%% Driftfusion 2D
ns1_df = u(:, 200, 2);
ps1_df = u(:, 261, 3);
ns2_df = u(:, 660, 2);
ps2_df = u(:, 721, 3);

R1_df2D = (ns1_df.*ps1_df - ni^2)./((1/sn1).*(ps1_df + pt1) + (1/sp1).*(ns1_df + nt1));
R2_df2D = (ns2_df.*ps2_df - ni^2)./((1/sn2).*(ps2_df + pt2) + (1/sp2).*(ns2_df + nt2));

%% Driftfusion 3D
rx = dfana.calcr(sol_df);

R1_df3D = zeros(length(t), 1);
for i=1:length(t)
    R1_df3D(i) = trapz(x(200:261), rx.srh(i, 200:261), 2);
    R2_df3D(i) = trapz(x(660:721), rx.srh(i, 660:721), 2);
end

%% IonMonger
t_im = sol_im.time;
ns1_im = sol_im.dstrbns.nE(:, end)*1e-6;
ps1_im = sol_im.dstrbns.p(:, 1)*1e-6;
ns2_im = sol_im.dstrbns.n(:, end)*1e-6;
ps2_im = sol_im.dstrbns.pH(:, 1)*1e-6;

R1_im = (ns1_im.*ps1_im - ni^2)./((1/sn1).*(ps1_im + pt1) + (1/sp1).*(ns1_im + nt1));
R2_im = (ns2_im.*ps2_im - ni^2)./((1/sn2).*(ps2_im + pt2) + (1/sp2).*(ns2_im + nt2));

%% Plots
figure(200)
semilogy(t, R1_df2D, t, R1_df3D, t_im, R1_im)
xlabel('Time [s]')
ylabel('Recombination flux [cm-2s-1]')
legend('DF 2D', 'DF volumetric', 'IM 2D')

figure(201)
semilogy(t, R2_df2D, t, R2_df3D, t_im, R2_im)
xlabel('Time [s]')
ylabel('Recombination flux [cm-2s-1]')
legend('DF 2D', 'DF volumetric', 'IM 2D')

end
