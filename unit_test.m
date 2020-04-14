%UNIT_TEST - Runs various simulations in order to test the functioning of Driftfusion and many other implemented functions
% For profiling the time spent running each function run 'profile on' before running the tests and
% 'profile viewer' after.
% Using Matlab's Coverage Reports, the obsolete and unused code can be easily spotted.
%
% Syntax:  runtests('unit_test')
%
% Other m-files required: pc, equilibrate, doJV, explore
% Subfunctions: none
% MAT-files required: none
%
% See also df.
%
% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

% prepare solutions and indirectly test equilibrate and genIntStructs
% functions
initialise_df
par = pc('Input_files/ptpd_mapi_pcbm.csv');
soleq = equilibrate(par);
JVsol = doJV(soleq.ion, 1e-2, 50, 1, 1, 0, 1.2, 3);
%[Ecb, Evb, Efn, Efp] = dfana.QFLs(soleq.ion);

% example taken from Scripts/explore_script
exsol = explore.explore2par(par, {'d(1,3)','Int'},...
    {[400e-7, 800e-7], logspace(-2,0,2)}, 20);

%% Core df no input
% solstruct = df(varargin)

df();
 
%% Core df one input
% solstruct = df(varargin)

df(soleq.ion);

%% Core df two inputs zero u
% solstruct = df(varargin)

sol = struct();
sol.u = 0;
df(sol, par);

%% Core df two inputs char
% solstruct = df(varargin)

df(soleq.ion, 'test');

%% Core df two inputs else
% solstruct = df(varargin)

df(soleq.ion, par);

%% Core dfplot ELx
% ELx(varargin)

dfplot.ELx(soleq.ion)
dfplot.ELx(JVsol.ill.f,[0,100])

%% Core dfplot Jt
% Jt(sol, xpos)

dfplot.Jt(JVsol.ill.f, 2e-5)

%% Core dfplot Jx
% Jx(varargin)

dfplot.Jx(soleq.ion)
dfplot.Jx(JVsol.ill.f,[0,100])

%% Core dfplot jx
% jx(varargin)

dfplot.jx(soleq.ion)
dfplot.jx(JVsol.ill.f,[0,100])

%% Core dfplot JV
% JV(JV, option)

dfplot.JV(JVsol, 3)

%% Core dfplot Jddx
% Jddx(varargin)

dfplot.Jddx(soleq.ion)
dfplot.Jddx(JVsol.ill.f,[0,100])

%% Core dfplot Voct
% Voct(sol)

dfplot.Voct(JVsol.ill.f)

%% Core dfplot PLt
% PLt(sol)

dfplot.PLt(JVsol.ill.f)

%% Core dfplot Vappt
% Vappt(sol)

dfplot.Vappt(JVsol.ill.f)

%% Core dfplot JVapp
% JVapp(sol, xpos)

dfplot.JVapp(soleq.ion, 2e-5)

%% Core dfplot JtotVapp
% JtotVapp(sol, xpos)

dfplot.JtotVapp(JVsol.ill.f, 2e-5)

%% Core dfplot logJVapp
% logJVapp(sol, xpos)

dfplot.logJVapp(JVsol.ill.f, 2e-5)

%% Core dfplot logJVapp3D
% logJVapp3D(sol, xpos, ylogon)
dfplot.logJVapp3D(JVsol.ill.f, 2e-5, true)

%% Core dfplot xmesh
% xmesh(sol)

dfplot.xmesh(soleq.ion)

%% Core dfplot Vx
% Vx(varargin)

dfplot.Vx(soleq.ion)
dfplot.Vx(JVsol.ill.f,[0,100])

%% Core dfplot npx
% npx(varargin)

dfplot.npx(soleq.ion)
dfplot.npx(JVsol.ill.f,[0,100])

%% Core dfplot acx
% acx(varargin)

dfplot.acx(soleq.ion)
dfplot.acx(JVsol.ill.f,[0,100])

%% Core dfplot gx
% gx(varargin)

dfplot.gx(JVsol.ill.f)
dfplot.gx(JVsol.ill.f,[0,100])

%% Core dfplot gxt
% gxt(sol)

dfplot.gxt(JVsol.ill.f)

%% Core dfplot rx
% rx(varargin)

dfplot.rx(JVsol.ill.f)
dfplot.rx(JVsol.ill.f,[0,100])

%% Core dfplot JVrec
% JVrec(JV, option)

dfplot.JVrec(JVsol, 3)

%% Core dfplot Ft
% Ft(sol, xpos)

dfplot.Ft(JVsol.ill.f, 2e-5)

%% Core dfplot sigmat
% sigmat(sol)

dfplot.sigmat(JVsol.ill.f)

%% Core dfplot Qt
% Qt(sol, x1, x2)

dfplot.Qt(JVsol.ill.f)

%% Core dfplot QVapp
% QVapp(sol, x1, x2)

dfplot.QVapp(JVsol.ill.f)

%% Core dfplot rhox
% rhox(varargin)

dfplot.rhox(soleq.ion)
dfplot.rhox(JVsol.ill.f,[0,100])

%% Core dfplot deltarhox
% deltarhox(varargin)

dfplot.deltarhox(JVsol.ill.f)
dfplot.deltarhox(JVsol.ill.f,[0,100])

%% Core dfplot rhoxFxVx
% rhoxFxVx(varargin)

dfplot.rhoxFxVx(soleq.ion)
dfplot.rhoxFxVx(JVsol.ill.f,[0,100])

%% Core dfplot rhoxVx
% rhoxVx(varargin)

dfplot.rhoxVx(soleq.ion)
dfplot.rhoxVx(JVsol.ill.f,[0,100])

%% Core dfplot ELx_single
% ELx_single(varargin)

dfplot.ELx_single(soleq.ion)
dfplot.ELx_single(JVsol.ill.f,[0,100])

%% Core dfplot ELnpx
% ELnpx(varargin)

dfplot.ELnpx(soleq.ion)
dfplot.ELnpx(JVsol.ill.f,[0,100])

%% Core dfplot Vacx
% Vacx(varargin)

dfplot.Vacx(soleq.ion)
dfplot.Vacx(JVsol.ill.f,[0,100])

%% Core dfplot Vionacx
% Vionacx(varargin)

dfplot.Vionacx(soleq.ion)
dfplot.Vionacx(JVsol.ill.f,[0,100])

%% Core dfplot Fiont
% Fiont(sol, xpos)

dfplot.Fiont(JVsol.ill.f, 2e-5)

%% Core dfplot PLx
% PLx(varargin)

dfplot.PLx(JVsol.ill.f)
dfplot.PLx(JVsol.ill.f,[0,100])

%% Core dfplot colourblocks
% colourblocks(sol, yrange)

dfplot.colourblocks(soleq.ion, [0,1])

%% Core dfplot sortarg
% [sol, tarr, pointtype, xrange] = sortarg(args)

test_tarr = [0,10,100];
test_xrange = [150, 250];
[~, ~, ~, ~] = dfplot.sortarg({JVsol.ill.f});
[~, tarr, ~, ~] = dfplot.sortarg({JVsol.ill.f, test_tarr});
assert(all(tarr == test_tarr))
[~, tarr, ~, xrange] = dfplot.sortarg({JVsol.ill.f, test_tarr, test_xrange});
assert(all([tarr == test_tarr, xrange == test_xrange]))

%% Core dfplot x2d
% x2d(sol, xmesh, variables, legstr, linestyle, ylab, tarr, xrange, logx, logy)

dfplot.x2d(JVsol.ill.f, JVsol.ill.f.x, {JVsol.ill.f.u(:,:,1)}, {'test'}, ['-','.'], 'ylab', JVsol.ill.f.t, [150,250], false, false)

%% Core build_device
% dev = build_device(par, meshoption)

dev1 = build_device(par, 'iwhole');
dev2 = build_device(par, 'ihalf');

%% Core build_property
% devprop = build_property(property, xmesh, par, interface_switch, gradient_property)

dev.taun = build_property(par.taun, par.xx, par, 'constant', 0);
dev.EA = build_property(par.EA, par.xx, par, 'lin_graded', 0);
dev.NA = build_property(par.NA, par.xx, par, 'log_graded', 0);
dev.g0 = build_property(par.g0, par.xx, par, 'zeroed', 0);
dev.gradEA = build_property(par.EA, par.xx, par, 'lin_graded', 1);
dev.gradNc = build_property(par.Nc, par.xx, par, 'log_graded', 1);

%% Core dfana splitsol
% [u,t,x,par,dev,n,p,a,c,V] = splitsol(sol)

[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(soleq.ion);

%% Core dfana QFLs
% [Ecb, Evb, Efn, Efp] = QFLs(sol)

[Ecb, Evb, Efn, Efp] = dfana.QFLs(soleq.ion);

%% Core dfana QFL_ihalf
% [Ecb, Evb, Efn, Efp] = QFL_ihalf(sol)

[Ecb, Evb, Efn, Efp] = dfana.QFL_ihalf(soleq.ion);

%% Core dfana QFL_J
% [Ecb, Evb, Efn, Efp] = QFL_J(sol)

[Ecb, Evb, Efn, Efp] = dfana.QFL_J(soleq.ion);

%% Core dfana calcJ
% [J, j, x] = calcJ(sol)

[J, j, x] = dfana.calcJ(soleq.ion);

%% Core dfana calcg
% [g1, g2, g] = calcg(sol)

[g1, g2, g] = dfana.calcg(soleq.ion);

%% Core dfana calcr
% r = calcr(sol)

r = dfana.calcr(soleq.ion);

%% Core dfana calcr_ihalf
% r = calcr_ihalf(sol)

r = dfana.calcr_ihalf(soleq.ion);

%% Core dfana Jddxt
% [jdd, Jdd, xout] = Jddxt(sol)

[jdd, Jdd, xout] = dfana.Jddxt(soleq.ion);

%% Core dfana calcF
% [FV, Frho] = calcF(sol)

[FV, Frho] = dfana.calcF(soleq.ion);

%% Core dfana calcF_ihalf
% [FV, Frho] = calcF_ihalf(sol)

[FV, Frho] = dfana.calcF_ihalf(soleq.ion);

%% Core dfana calcrho
% rho = calcrho(sol)

rho = dfana.calcrho(soleq.ion);

%% Core dfana calcrho_ihalf
% rho = calcrho_ihalf(sol)

rho2 = dfana.calcrho_ihalf(soleq.ion);

%% Core dfana calcVapp
% Vapp = calcVapp(sol)

Vapp = dfana.calcVapp(soleq.ion);

%% Core dfana JVstats
% stats = JVstats(JVsol)

stats = dfana.JVstats(JVsol);

%% Core dfana PLt
% value = PLt(sol)

value = dfana.PLt(soleq.ion);

%% Core dfana calcVQFL
% VQFL = calcVQFL(sol)

VQFL = dfana.calcVQFL(soleq.ion);

%% Core dfana deltaVt
% deltaV = deltaVt(sol, p1, p2)

deltaV = dfana.deltaVt(soleq.ion, 1, 123);

%% Core dfana calcsigma
% sigma = calcsigma(sol)

sigma = dfana.calcsigma(soleq.ion);

%% Core dfana calcsigma_ion
% sigma_ion = calcsigma_ion(sol)

sigma_ion = dfana.calcsigma_ion(soleq.ion);

%% Core dfana calcFion
% Fion = calcFion(sol)

Fion = dfana.calcFion(soleq.ion);

%% Core dfana calcVion
% Vion = calcVion(sol)

Vion = dfana.calcVion(soleq.ion);

%% Core dfana pdentrp
% [U,Ux] = pdentrp(singular,m,xL,uL,xR,uR,xout)

[Vloc,dVdx] = dfana.pdentrp(false,false,par.xx(123),V(1,123),par.xx(124),V(1,124),1e-6);

%% Core ditro_fun nfun
% n = nfun(Nc, Ec, Efn, T, prob_distro_function)

n = distro_fun.nfun(par.dev.Nc, par.dev.EA, par.dev.Et, par.T, par.prob_distro_function);

%% Core ditro_fun pfun
% p = pfun(Nv, Ev, Efp, T, prob_distro_function)

p = distro_fun.pfun(par.dev.Nv, par.dev.IP, par.dev.Et, par.T, par.prob_distro_function);

%% Core ditro_fun Dn_fd_fun and Dnlook and Efn_fd_fun
% Dnfd = Dn_fd_fun(Nc, Ec, Efn, mue, T)

Ec = -4.95;
Efn = -5;
Dnfd = distro_fun.Dn_fd_fun(par.dev.Nc(end), Ec, Efn, par.mue(end), par.T);

% Dsol = Dnlook(n, Dnfun, n_fd)

Dsol = distro_fun.Dnlook(soleq.ion.u(2,1,end), Dnfd.Dnfun, Dnfd.n_fd);

% Efn_fd = Efn_fd_fun(n, Efn, n_fd)

Efn_fd = distro_fun.Efn_fd_fun(soleq.ion.u(2,1,end), Efn, Dnfd.n_fd);

%% Core ditro_fun Dp_fd_fun and Dplook and Efp_fd_fun
% Dpfd = Dp_fd_fun(Nv, Ev, Efp, muh, T)

Dpfd = distro_fun.Dp_fd_fun(par.dev.Nv(1), par.IP(1), Efp, par.muh(1), par.T);

% Dsol = Dplook(p, Dpfun, p_fd)

Dsol = distro_fun.Dplook(soleq.ion.u(3,1,1), Dpfd.Dpfun, Dpfd.p_fd);

% Efp_fd = Efp_fd_fun(p, Efp, p_fd)

Efp_fd = distro_fun.Efp_fd_fun(soleq.ion.u(3,1,1), Efp, Dpfd.p_fd);

%% Core fun_gen
% fun = fun_gen(fun_type)

fun = fun_gen(fun_type);

%% Core generation
% gx = generation(par, source_type, laserlambda)

gx = generation(par, source_type, laserlambda);

%% Core getdevihalf
% devihalf = getdevihalf(par)

devihalf = getdevihalf(par);

%% Core getvarihalf
% varihalf = getvarihalf(var)

varihalf = getvarihalf(var);

%% Core getxihalf
% xsolver = getxihalf(sol)

xsolver = getxihalf(sol);

%% Core import_properties
% par = import_properties(par, filepath)

par = import_properties(par, filepath);

%% Core meshgen_t default
% [t] = meshgen_t(par)

part = par;
part.mesht_figon = true;
[t] = meshgen_t(part);

%% Core meshgen_t 1
% [t] = meshgen_t(par)

part = par;
part.mesht_figon = true;
part.tmesh_type = 1;
[t] = meshgen_t(part);

%% Core meshgen_t 2
% [t] = meshgen_t(par)

part = par;
part.mesht_figon = true;
part.tmesh_type = 2;
[t] = meshgen_t(part);

%% Core meshgen_t 3
% [t] = meshgen_t(par)

part = par;
part.mesht_figon = true;
part.tmesh_type = 3;
[t] = meshgen_t(part);

%% Core meshgen_t 4
% [t] = meshgen_t(par)

part = par;
part.mesht_figon = true;
part.tmesh_type = 4;
[t] = meshgen_t(part);

%% Core meshgen_x default
% x = meshgen_x(par)

parx = par;
parx.meshx_figon = true;
x0 = meshgen_x(par);

%% Core meshgen_x 1
% x = meshgen_x(par)

parx = par;
parx.meshx_figon = true;
parx.xmesh_type = 1;
x1 = meshgen_x(parx);

%% Core meshgen_x 2
% x = meshgen_x(par)

parx = par;
parx.meshx_figon = true;
parx.xmesh_type = 2;
x2 = meshgen_x(parx);

%% Core meshgen_x 3
% x = meshgen_x(par)

parx = par;
parx.meshx_figon = true;
parx.xmesh_type = 3;
x3 = meshgen_x(parx);

%% Core meshgen_x 4
% x = meshgen_x(par)

parx = par;
parx.meshx_figon = true;
parx.xmesh_type = 4;
x4 = meshgen_x(parx);

%% Core meshgen_x 5
% x = meshgen_x(par)

parx = par;
parx.meshx_figon = true;
parx.xmesh_type = 5;
x5 = meshgen_x(parx);

%% Core refresh_device
% par = refresh_device(par)

par = refresh_device(par);

%% Core triangle_fun
% y = triangle_fun(coeff, t)

y = triangle_fun(coeff, t);

%% Input_files

inputs = dir('Input_files');
for i=1:length(inputs)
    input = inputs(i).name;
    if ~any(regexp(input,'^\.'))
        par = pc(input);
        xpoints = round(1 + sum(par.layer_points));

        soleq = equilibrate(par);
        
        el = soleq.el;
        el_s = size(el.u);

        exp_el_s = [el.par.tpoints, xpoints, 3];
        assert(all(el_s == exp_el_s), 'Expected size: %d, %d, %d. Obtained size: %d, %d, %d.', el_s(1), el_s(2), el_s(3), exp_el_s(1), exp_el_s(2), exp_el_s(3));
        assert(~any(isnan(el.u(:))))
        
        ion = soleq.ion;
        ion_s = size(ion.u);
        
        exp_ion_s = [ion.par.tpoints, xpoints, round(3+par.N_ionic_species)];
        assert(all(ion_s == exp_ion_s), 'Expected size: %d, %d, %d. Obtained size: %d, %d, %d.', el_s(1), el_s(2), el_s(3), exp_el_s(1), exp_el_s(2), exp_el_s(3));
        assert(~any(isnan(ion.u(:))))
    end
end

%% Scripts

inputs = dir('Scripts');
for i=1:length(inputs)
    input = inputs(i).name;
    if ~any(regexp(input,'^\.|^test'))
        
        disp(['### running script ' input]);
        run(input);
    end
end

%% Helper calcJsc, calcR0 and Eg_vs_Voc
% [EgArr, Jsc_vs_Eg] = calcJsc
[EgArr, Jsc_vs_Eg] = calcJsc;

% [JV_ana, r0, k_rad, Voc, g0] = calcR0(EgArr, Jsc_vs_Eg, par)
[JV_ana, r0, k_rad, Voc, g0] = calcR0(EgArr, Jsc_vs_Eg, par);

% [G0_Arr, k_rad_Arr, R0_Arr, Eg, VocArr, DeltaVoc] = Eg_vs_Voc(EgArr, Jsc_vs_Eg)
[G0_Arr, k_rad_Arr, R0_Arr, Eg, VocArr, DeltaVoc] = Eg_vs_Voc(EgArr, Jsc_vs_Eg);

%% Helper explore plotPL
% plotPL(exsol)
explore.plotPL(exsol);

%% Helper explore plotsurf
% plotsurf(exsol, yproperty, xlogon, ylogon, zlogon)
explore.plotsurf(exsol, 'Voc_r', true, false, false)

%% Helper explore getJtot
% Jtot = getJtot(sol)
Jtot = explore.getJtot(soleq.ion);

%% Helper explore writevar
% var = writevar(var, i, j, xx, arr)
x=0;
x = explore.writevar(x, 1, 1, par.xx, soleq.ion.x);

%% Helper explore helper
% par = helper(par, parname, parvalue)
par2 = explore.helper(par, 'Rs', 1);

%% Helper explore plotstat_2D_parval1
% plotstat_2D_parval1(exsol, yproperty, logx, logy)
explore.plotstat_2D_parval1(exsol, 'Voc_r', true, false)

%% Helper explore plotstat_2D_parval2
% plotstat_2D_parval2(exsol, yproperty, logx, logy)
explore.plotstat_2D_parval2(exsol, 'Voc_r', true, false)

%% Helper explore plotfinalELx
% plotfinalELx(exsol)
explore.plotfinalELx(exsol)

%% Helper explore plotprof_2D
% plotprof_2D(exsol, yproperty, par1logical, par2logical, logx,logy)
explore.plotprof_2D(exsol, 'Voc_r', true, false)

%% Helper explore plotU
% plotU(exsol, par1logical, par2logical,logx,logy)
explore.plotU(exsol, [true,true], [true,true],false,false)

%% Helper explore plotCE
% plotCE(exsol_Voc, exsol_eq, xlogon, ylogon, zlogon, normalise)
explore.plotCE(exsol, exsol, false, false, false, "ciaomamma")

%% Helper explore plotJV
% plotJV(exsol, par1logical, par2logical)
explore.plotJV(exsol, [true,true], [true,true])

%% Helper getpointpos
% ppos = getpointpos(xpos, xmesh)
ppos = getpointpos(1e-6, soleq.ion.x);

%% Helper importsolcoresol
% solstruct = importsolcoresol(varargin)
solstruct = importsolcoresol();

%% Helper makemovie
% Framefile = makemovie(sol, plotfun, xrange, yrange, movie_name, Vcounter, tcounter)
Framefile = makemovie(JVsol.ill.f, @dfplot.PLx, [0,100e-7], 0, 'test_makemovie_delete_me', true, true);

%% Helper verifyStabilization
% all_stable = verifyStabilization(sol_matrix, t_array, time_fraction)
all_stable = verifyStabilization(soleq.ion.u, soleq.ion.t, 0.1);

%% Protocols changeLight
% sol_int = changeLight(sol, newInt, tmax)
sol_int1 = changeLight(soleq.ion, 1, 2);

%% Protocols doCV
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV = doCV(soleq.ion, 0.1, 0, 0.5, -0.5, 0.1, 1, 50);

%% Protocols doIMPS
% sol_IMPS = doIMPS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
sol_IMPS = doIMPS(soleq.ion, 0.2, 0.2, 1, 10, 50);

%% Protocols doIMVS
% sol_IMVS = doIMVS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
sol_IMVS = doIMVS(soleq.ion, 0.2, 0.2, 1, 10, 50);

%% Protocols doLightPulse
% [sol_pulse] = doLightPulse(sol_ini, pulse_int, tmax, tpoints, duty, mobseti, log_timemesh)
sol_pulse = doLightPulse(soleq.ion, 0.1, 1, 20, 5, true, true);
sol_pulse2 = doLightPulse(soleq.ion, 0.1, 1, 20, 5, true, false);

%% Protocols doSDP and Analysis anasdp
% sdpsol = doSDP(sol_ini, tdwell_arr, Vjump, pulse_int, pulse_tmax, duty, tpoints, scalefactor)
sdpsol = doSDP(soleq.ion, [1e-3,1], 0.9, 1, 1, 5, 50, 1);

% anasdp(sdpsol, Jtr_time)
anasdp(sdpsol, 1);

%% Protocols doSPV and Analysis spvana
% spvsol = doSPV(sol_ini, Int, mobseti, tpoints, tmax, Rs, stabilise)
spvsol = doSPV(soleq.ion, 1, true, 20, 1, 1e3, true);

% spvdat = spvana(spvsol)
spvdat = spvana(spvsol);

%% Protocols doTPV
% [sol_TPV, sol_ill] = doTPV(sol_ini, bias_int, stab_time, mobseti, Rs, pulse_int, tmax, tpoints, duty)
[sol_TPV, sol_ill] = doTPV(soleq.ion, 1, 10, true, 1e-2, 0.1, 1, 20, 5);

%% Protocols findVocDirect
% [sol_Voc, Voc] = findVocDirect(sol_ini, light_intensity, mobseti)
[sol_Voc, Voc] = findVocDirect(soleq.ion, 1, true);

%% Protocols findVoc
% [sol_Voc, Voc] = findVoc(sol_ini, Int, mobseti, x0, x1, tol)
[sol_Voc, Voc] = findVoc(soleq.ion, 1, true, 0.9, 1.3, 1e-5)

%% Protocols genIntStructs
% [structCell, V_array, J_array] = genIntStructs(struct_eq, startInt, endInt, points, include_dark)
[structCell, V_array, J_array] = genIntStructs(soleq.ion, 1e-4, 1, 5, true);

%% Protocols genIntStructsRealVoc
% [goodVocAsymStructCell, VOCs] = genIntStructsRealVoc(struct_eq, startInt, endInt, points, include_dark)
[goodVocAsymStructCell, VOCs] = genIntStructsRealVoc(soleq.ion, 1e-4, 1e-3, 2, false);

%% Protocols genVappStructs
% VappSol = genVappStructs(solini, Vapp_arr, mobseti)
VappSol = genVappStructs(soleq.ion, [0.9,1.1,1.3], true);

%% Protocols jumptoV
% sol_relax = jumptoV(sol_ini, Vjump, tdwell, mobseti, Int, stabilise, accelerate)
sol_relax = jumptoV(soleq.ion, 1.2, 1, true, 1, true, true);

%% Protocols lightonRs
% sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
sol_ill = lightonRs(soleq.ion, 1, 10, true, 1e6, 50);

%% Protocols stabilize
% steadystate_struct = stabilize(struct)
steadystate_struct = stabilize(struct);

%% Protocols sweeplight
% sol_sweep = sweeplight(sol_ini, tmax, tpoints, end_int)
sol_sweep = sweeplight(soleq.ion, 1, 50, 1);

%% Protocols transient_nid and Analysis transient_nid_ana
% sol_OC = transient_nid(sol_ini, int_arr, stab_time, mobseti, Rs, pnts)
sol_OC = transient_nid(soleq.ion, [0.01,0.1,1], 10, true, 1e9, 10);
nidt = transient_nid_ana(sol_OC);

%% Protocols VappFunction
% sol = VappFunction(sol_ini, Vapp_func, Vapp_coeff, tmax, tpoints, logtime)
sol = VappFunction(soleq.ion, 'constant', 0.2, 10, 30, false);
sol2 = VappFunction(soleq.ion, 'sweep', [0,1,5], 10, 30, false);
sol3 = VappFunction(soleq.ion, 'square', [0,1,5,30], 10, 30, false);
sol4 = VappFunction(soleq.ion, 'sin', [0, 1e-2, 1, 0], 10, 30, false);

%% Optical beerlambert
% Gentot = beerlambert(par, x, source_type, laserlambda, figson)
Gentot = beerlambert(par, par.xx, 'AM15', 0, true);
Gentot2 = beerlambert(par, par.xx, 'laser', 500, true);

%% Optical lightsource
% AM15 = lightsource(source_type, lambda)
AM15 = lightsource('AM15', 500);

%% Optical LoadRefrIndex
% [n_interp, k_interp] = LoadRefrIndex(name,wavelengths)
[n_interp, k_interp] = LoadRefrIndex(par.stack{1},300:767);

%------------- END OF CODE --------------

