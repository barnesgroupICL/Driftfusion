function compare_rec_flux(sol_df, sol_im, plot_switch, tolerance)
% Script to compare the interfacial recombination fluxes from Driftfusion
% (DF) and IonMonger (IM)
% Currently only working for ETL-AL-HTL architecture
% SOL_IM is either a solution from IonMonger or if SOL_IM = 0, then only
% the DF vectors are plotted.
% PLOT_SWITCH = 1 for plot outputs
% TOLERANCE is the fractional difference between DF 2D and volumetric calculations above
% which a warning is displayed
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

%% Driftfusion Abrupt
ns1_df = u(:, par.pcum0(2)+1, 2);
ps1_df = u(:, par.pcum0(3)+1, 3);
ns2_df = u(:, par.pcum0(4)+1, 2);
ps2_df = u(:, par.pcum0(5)+1, 3);

R1_df2D = (ns1_df.*ps1_df - ni^2)./((1/sn1).*(ps1_df + pt1) + (1/sp1).*(ns1_df + nt1));
R2_df2D = (ns2_df.*ps2_df - ni^2)./((1/sn2).*(ps2_df + pt2) + (1/sp2).*(ns2_df + nt2));

%% Driftfusion 3D
rx = dfana.calcr_ihalf(sol_df);

R1_df3D = zeros(length(t), 1);
R2_df3D = zeros(length(t), 1);
for i=1:length(t)
    R1_df3D(i) = trapz(x(par.pcum0(2)+1:par.pcum0(3)+1), rx.vsr(i, par.pcum0(2)+1:par.pcum0(3)+1), 2);
    R2_df3D(i) = trapz(x(par.pcum0(4)+1:par.pcum0(5)+1), rx.vsr(i, par.pcum0(4)+1:par.pcum0(5)+1), 2);
end

%% Fractional difference
delta1 = 1-(R1_df3D./R1_df2D);
delta2 = 1-(R2_df3D./R2_df2D);

if [max(delta1), max(delta2)] > tolerance
    warning(['The volumetric surface recombination model error exceeded the user-defined tolerance level (tol_vsr = ', ...
        num2str(tolerance), '). Consider increasing the interface layer electronic mobilities or switching to a alternative recombination model.'])
end

%% IonMonger
if sol_im ~= 0
    t_im = sol_im.time;
    ns1_im = sol_im.dstrbns.nE(:, end)*1e-6;
    ps1_im = sol_im.dstrbns.p(:, 1)*1e-6;
    ns2_im = sol_im.dstrbns.n(:, end)*1e-6;
    ps2_im = sol_im.dstrbns.pH(:, 1)*1e-6;
    
    R1_im = (ns1_im.*ps1_im - ni^2/sol_im.params.kE)./((1/sn1).*(ps1_im + pt1) + (1/sp1).*(ns1_im + nt1));
    R2_im = (ns2_im.*ps2_im - ni^2/sol_im.params.kH)./((1/sn2).*(ps2_im + pt2) + (1/sp2).*(ns2_im + nt2));
end

if plot_switch
    figure(200)
    semilogy(t, R1_df2D, t, R1_df3D)
    if sol_im ~= 0
        hold on
        semilogy(t_im, R1_im)
        hold off
    end
    xlabel('Time [s]')
    ylabel('Recombination flux [cm-2s-1]')
    legend('DF 2D', 'DF volumetric', 'IM 2D')
    
    figure(201)
    semilogy(t, R2_df2D, t, R2_df3D)
    if sol_im ~= 0
        hold on
        semilogy(t_im, R2_im)
        hold off
    end
    xlabel('Time [s]')
    ylabel('Recombination flux [cm-2s-1]')
    legend('DF 2D', 'DF volumetric', 'IM 2D')
    
    figure(202)
    semilogy(t, ns1_df, t, ps1_df)
    if sol_im ~= 0
        hold on
        semilogy(t_im, ns1_im, '--', t_im, ps1_im, '--')
        hold off
    end
    xlabel('Time [s]')
    ylabel('Carrier density')
    legend('ns1 DF', 'ps1 DF', 'ns1 IM', 'ps1 IM')
    
    figure(203)
    semilogy(t, ns2_df, t, ps2_df)
    if sol_im ~= 0
        hold on
        semilogy(t_im, ns2_im, '--', t_im, ps2_im, '--')
        hold off
    end
    xlabel('Time [s]')
    ylabel('Carrier density')
    legend('ns2 DF', 'ps2 DF', 'ns2 IM', 'ps2 IM')
    
    figure(204)
    plot(t, delta1, t, delta2)
    xlabel('Time [s]')
    ylabel('Fractional difference')
    legend('interface 1', 'interface 2')
end
end
