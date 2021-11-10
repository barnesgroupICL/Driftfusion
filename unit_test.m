% UNIT_TEST - Runs various simulations in order to test the functioning of Driftfusion and many other implemented functions
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
% clear all
% prepare solutions and indirectly test equilibrate and genIntStructs
% functions
initialise_df

% for testing more realistic parameters
%input_csv = 'Input_files/pedotpss_mapi_pcbm.csv';
% a structure written for testing purposes
input_csv = 'Input_files/3_layer_unit_test.csv';

% par = pc(varargin)
par = pc(input_csv);
% Relax VSR tolerance to 10% for test case to avoid warnings
par.RelTol_vsr = 0.1;
% soleq = equilibrate(varargin)
soleq = equilibrate(par);
% JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
JVsol = doJV(soleq.ion, 1e-2, 50, 1, true, 0, 1.0, 3);

%% Core df no input

% solstruct = df(varargin)
df();
 
%% Core df one input

% solstruct = df(varargin)
df(soleq.ion);

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
dfplot.JrecVapp(JVsol, 3)

%% Core dfplot Ft

% Ft(sol, xpos)
dfplot.Ft(JVsol.ill.f, 2e-5)

%% Core dfplot sigmat

% sigmat(sol)
dfplot.sigmat(JVsol.ill.f)

%% Core dfplot Qt

% Qt(sol, x1, x2)
dfplot.Qt(JVsol.ill.f, JVsol.ill.f.x(1), JVsol.ill.f.x(end))

%% Core dfplot QVapp

% QVapp(sol, x1, x2)
dfplot.QVapp(JVsol.ill.f, JVsol.ill.f.x(1), JVsol.ill.f.x(end))

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
dfplot.ELx(soleq.ion)
dfplot.ELx(JVsol.ill.f,[0,100])

%% Core dfplot ELnpx

% ELnpx(varargin)
dfplot.ELnpx(soleq.ion)
dfplot.ELnpx(JVsol.ill.f,[0,100])

%% Core dfplot Vxacx

% Vacx(varargin)
dfplot.Vxacx(soleq.ion)
dfplot.Vxacx(JVsol.ill.f,[0,100])

%% Core dfplot Vionacx

% Vionacx(varargin)
dfplot.Vionxacx(soleq.ion)
dfplot.Vionxacx(JVsol.ill.f,[0,100])

%% Core dfplot Fiont

% Fiont(sol, xpos)
dfplot.Fiont(JVsol.ill.f, 2e-5)

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
dev1 = build_device(par, 'whole');
dev2 = build_device(par, 'sub');

%% Core build_property

% devprop = build_property(property, xmesh, par, interface_switch, gradient_property)
build_property(par.taun, par.xx, par, 'constant', 0);
build_property(par.Phi_EA, par.xx, par, 'lin_graded', 0);
build_property(par.NA, par.xx, par, 'log_graded', 0);
build_property(par.g0, par.xx, par, 'zeroed', 0);
build_property(par.Phi_EA, par.xx, par, 'lin_graded', 1);
build_property(par.Nc, par.xx, par, 'log_graded', 1);

%% Core dfana splitsol

% [u,t,x,par,dev,n,p,a,c,V] = splitsol(sol)
dfana.splitsol(soleq.ion);

%% Core dfana calcEnergies

% [Ecb, Evb, Efn, Efp] = calcEnergies(sol)
dfana.calcEnergies(soleq.ion);

%% Core dfana calcJ

% [J, j, x] = calcJ(sol)
dfana.calcJ(soleq.ion);

%% Core dfana calcg

% [g1, g2, g] = calcg(sol)
dfana.calcg(soleq.ion);

%% Core dfana calcr

% r = calcr(sol, mesh_option)
dfana.calcr(soleq.ion, "whole");
dfana.calcr(soleq.ion, "sub");

%% Core dfana calcJdd

% [Jdd, jdd, xout] = calcJdd(sol)
dfana.calcJdd(soleq.ion);

%% Core dfana calcF

% [FV, Frho] = calcF(sol, mesh_option)
dfana.calcF(soleq.ion, "whole");
dfana.calcF(soleq.ion, "sub");

%% Core dfana calcrho

% rho = calcrho(sol, mesh_option)
dfana.calcrho(soleq.ion, "whole");
dfana.calcrho(soleq.ion, "sub");

%% Core dfana calcVapp

% Vapp = calcVapp(sol)
dfana.calcVapp(soleq.ion);

%% Core dfana JVstats

% stats = JVstats(JVsol)
dfana.JVstats(JVsol);

%% Core dfana calcPLt

% value = calcPLt(sol)
dfana.calcPLt(soleq.ion);

%% Core dfana calcDeltaQFL

% VQFL = calcDeltaQFL(sol)
dfana.calcDeltaQFL(soleq.ion);

%% Core dfana deltaVt

% deltaV = deltaVt(sol, p1, p2)
dfana.deltaVt(soleq.ion, 1, 123);

%% Core dfana calcsigma

% sigma = calcsigma(sol)
dfana.calcsigma(soleq.ion);

%% Core dfana calcsigma_ion

% sigma_ion = calcsigma_ion(sol)
dfana.calcsigma_ion(soleq.ion);

%% Core dfana calcFion

% Fion = calcFion(sol)
dfana.calcFion(soleq.ion);

%% Core dfana calcVion

% Vion = calcVion(sol)
dfana.calcVion(soleq.ion);

%% Core dfana pdentrp

% [U,Ux] = pdentrp(singular,m,xL,uL,xR,uR,xout)
dfana.pdentrp(false,false,par.xx(123),soleq.ion.u(1,123,1),par.xx(124),soleq.ion.u(1,124,1),1e-6);

%% Core ditro_fun nfun

% n = nfun(Nc, Ec, Efn, T, prob_distro_function)
distro_fun.nfun(par.dev.Nc, par.dev.Phi_EA, par.dev.EF0, par);

%% Core ditro_fun pfun

% p = pfun(Nv, Ev, Efp, T, prob_distro_function)
distro_fun.pfun(par.dev.Nv, par.dev.Phi_IP, par.dev.EF0, par);

%% Core ditro_fun Dn_fd_fun and Dnlook and Efn_fd_fun

Ec = -4.95;
[~, ~, Efn, ~] = dfana.calcEnergies(soleq.ion);
% Dnfd = Dn_fd_fun(Nc, Ec, Efn, mu_n, T)
Dnfd = distro_fun.Dn_fd_fun(par.dev.Nc(end), Ec, Efn, par.mu_n(end), par.T);

% Dsol = Dnlook(n, Dnfun, n_fd)
distro_fun.Dnlook(soleq.ion.u(2,1,end), Dnfd.Dnfun, Dnfd.n_fd);

% Efn_fd = Efn_fd_fun(n, Efn, n_fd)
distro_fun.Efn_fd_fun(soleq.ion.u(2,1,end), Efn, Dnfd.n_fd);

%% Core ditro_fun Dp_fd_fun and Dplook and Efp_fd_fun

[~, ~, ~, Efp] = dfana.calcEnergies(soleq.ion);
% Dpfd = Dp_fd_fun(Nv, Ev, Efp, mu_p, T)
Dpfd = distro_fun.Dp_fd_fun(par.dev.Nv(1), par.Phi_IP(1), Efp, par.mu_p(1), par.T);

% Dsol = Dplook(p, Dpfun, p_fd)
distro_fun.Dplook(soleq.ion.u(3,1,1), Dpfd.Dpfun, Dpfd.p_fd);

% Efp_fd = Efp_fd_fun(p, Efp, p_fd)
distro_fun.Efp_fd_fun(soleq.ion.u(3,1,1), Efp, Dpfd.p_fd);

%% Core fun_gen

% fun = fun_gen(fun_type)
fun_gen('constant');
fun_gen('sweep');
fun_gen('square');
fun_gen('sin');
fun_gen('tri');

%% Core generation
% gx = generation(par, source_type, laserlambda)
par_om1 = par;
par_om1.optical_model = 'uniform';
generation(par_om1, 'AM15', 470);
generation(par_om1, 'laser', 470);

par_om2 = par;
par_om2.optical_model = 'Beer-Lambert';
generation(par_om2, 'AM15', 470);
generation(par_om2, 'laser', 470);

%% Core getvar_sub

% varsub = getvar_sub(var)
EA_sub = getvar_sub(par.dev.Phi_EA);

%% Core getx_sub

% xsolver = getx_sub(sol)
getx_sub(soleq.ion);

%% Core import_properties

% par = import_properties(par, filepath)
import_properties(par, {input_csv});

%% Core meshgen_t default

% [t] = meshgen_t(par)
part = par;
meshgen_t(part);

%% Core meshgen_t 1

% [t] = meshgen_t(par)
part = par;
part.tmesh_type = 'linear';
meshgen_t(part);

%% Core meshgen_t 2

% [t] = meshgen_t(par)
part = par;
part.tmesh_type = 'log10';
meshgen_t(part);

%% Core meshgen_t 3

% [t] = meshgen_t(par)
part = par;
part.tmesh_type = 'log10-double';
meshgen_t(part);

%% Core meshgen_x default

% x = meshgen_x(par)
parx = par;
meshgen_x(parx);

%% Core meshgen_x 4

% x = meshgen_x(par)
parx = par;
parx.xmesh_type = 'linear';
meshgen_x(parx);

%% Core meshgen_x 5

% x = meshgen_x(par)
parx = par;
parx.xmesh_type = 'erf-linear';
meshgen_x(parx);

%% Core refresh_device

% par = refresh_device(par)
refresh_device(par);

%% Core triangle_fun

% y = triangle_fun(coeff, t)
triangle_fun([0.1, 0.2, 1, 3, 5], 0:0.1:10);

%% Analysis
% sigma_sum_filter = compare_rec_flux(sol_df, RelTol_vsr, AbsTol_vsr, plot_switch)
sigma_sum_R_flux = compare_rec_flux(JVsol.ill.f, 1e6, 0.05, 1);

%% Analysis
% [n_ana, p_ana, jn_ana, jp_ana] = compare_carrier_interfaces(sol, tarr, plot_switch)
[n_ana, p_ana, jn_ana, jp_ana] = compare_carrier_interfaces(JVsol.ill.f, JVsol.ill.f.t(end)*[0, 0.2, 0.4, 0.6], 1);

%% Helper calcJsc, calcR0 and Eg_vs_Voc

% [EgArr, Jsc_vs_Eg] = calcJsc
[EgArr, Jsc_vs_Eg] = calcJsc;

% [JV_ana, r0, k_rad, Voc, g0] = calcR0(EgArr, Jsc_vs_Eg, par)
calcR0(EgArr, Jsc_vs_Eg, par);

%% Helper getpointpos

% ppos = getpointpos(xpos, xmesh)
getpointpos(1e-6, soleq.ion.x);

%% Helper makemovie

% Framefile = makemovie(sol, plotfun, xrange, yrange, movie_name, Vcounter, tcounter)
makemovie(JVsol.dk.f, @dfplot.ELx, [0,100e-7], 0, 'test_makemovie_delete_me', true, true);

%% Helper verifyStabilization

% all_stable = verifyStabilization(sol_matrix, t_array, time_fraction)
verifyStabilization(soleq.ion.u, soleq.ion.t, 0.1);

%% Protocols changeLight

% sol_int = changeLight(sol, newInt, tmax)
changeLight(soleq.ion, 1, 2);

%% Protocols doCV

% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
doCV(soleq.ion, 0.1, 0, 0.5, -0.5, 0.1, 1, 50);

%% Protocols doIMPS

% sol_IMPS = doIMPS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
doIMPS(soleq.ion, 0.2, 0.02, 1, 10, 50);

%% Protocols doIMVS

% sol_IMVS = doIMVS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
doIMVS(soleq.ion, 1, 0.2, 4, 1, 50);

%% Protocols doLightPulse

% [sol_pulse] = doLightPulse(sol_ini, pulse_int, tmax, tpoints, duty, mobseti, log_timemesh)
doLightPulse(soleq.ion, 0.1, 1, 20, 5, true, true);
doLightPulse(soleq.ion, 0.1, 1, 20, 5, true, false);

%% Protocols doSDP and Analysis anasdp

% sdpsol = doSDP(sol_ini, tdwell_arr, Vjump, pulse_int, pulse_tmax, duty, tpoints, scalefactor)
sdpsol = doSDP(soleq.ion, [1e-3,1], 0.9, 1, 1, 5, 50, 1);

% anasdp(sdpsol, Jtr_time)
anasdp(sdpsol, 1);

%% Protocols doSPV and Analysis spvana

% spvsol = doSPV(sol_ini, Int, mobseti, tpoints, tmax, Rs, stabilise)
spvsol = doSPV(soleq.ion, 1, true, 20, 1, 1e3, true);

% spvdat = spvana(spvsol)
spvana(spvsol);

%% Protocols doTPV

% [sol_TPV, sol_ill] = doTPV(sol_ini, bias_int, stab_time, mobseti, Rs, pulse_int, tmax, tpoints, duty)
doTPV(soleq.ion, 1, 10, true, 1e-2, 0.1, 1, 20, 5);

%% Protocols findVocDirect

% [sol_Voc, Voc] = findVocDirect(sol_ini, light_intensity, mobseti, tpoints)
findVocDirect(soleq.ion, 1, true, 40);

%% Protocols findVoc

% [sol_Voc, Voc] = findVoc(sol_ini, Int, mobseti, x0, x1, tol, tpoints, plot)
findVoc(soleq.ion, 1, true, 0.9, 1.3, 1e-5, 50, 0)

%% Protocols genIntStructs

% [structCell, V_array, J_array] = genIntStructs(struct_eq, startInt, endInt, points, include_dark)
genIntStructs(soleq.ion, 1e-4, 1, 5, true);

%% Protocols genIntStructsRealVoc

% [goodVocAsymStructCell, VOCs] = genIntStructsRealVoc(struct_eq, startInt, endInt, points, include_dark)
genIntStructsRealVoc(soleq.ion, 1e-4, 1e-3, 2, false);

%% Protocols genVappStructs

% VappSol = genVappStructs(solini, Vapp_arr, mobseti)
genVappStructs(soleq.ion, [0.9,1.1,1.3], true);

%% Protocols jumptoV

% sol_relax = jumptoV(sol_ini, Vjump, tdwell, mobseti, Int, stabilise, accelerate)
jumptoV(soleq.ion, 1.2, 1, true, 1, true, true);

%% Protocols lightonRs

% sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
lightonRs(soleq.ion, 1, 10, true, 1e6, 50);

%% Protocols stabilize

% steadystate_struct = stabilize(struct)
stabilize(soleq.ion);

%% Protocols sweeplight

% sol_sweep = sweeplight(sol_ini, tmax, tpoints, end_int)
sweeplight(soleq.ion, 1, 50, 1);

%% Protocols transient_nid and Analysis transient_nid_ana

% sol_OC = transient_nid(sol_ini, int_arr, stab_time, mobseti, Rs, pnts)
sol_OC = transient_nid(soleq.ion, [0.01,0.1,1], 10, true, 1e9, 10);
transient_nid_ana(sol_OC);

%% Protocols VappFunction

% sol = VappFunction(sol_ini, Vapp_func, Vapp_coeff, tmax, tpoints, logtime)

% coeff(1);
VappFunction(soleq.ion, 'constant', 0, 10, 30, false);

% coeff(1) + (coeff(2)-coeff(1))*t/coeff(3);
VappFunction(soleq.ion, 'sweep', [0, 1, 5], 10, 30, false);

% coeff(1) + (coeff(2)-coeff(1))*lt(mod(t,coeff(3))*1/coeff(3),coeff(4)/100);
% VappFunction(soleq.ion, 'square', [0.1, 0, 10, 90], 10, 30, false);

% coeff(1) + coeff(2)*(sin(2*pi*coeff(3)*t + coeff(4)));
VappFunction(soleq.ion, 'sin', [0, 0.1, 1, 0], 10, 30, false);

% COEFF = [OFFSET, V1, V2, periods, tperiod]
% triangle_fun(coeff, t);
VappFunction(soleq.ion, 'tri', [0, 0.6, -0.1, 3, 5], 15, 60, false);

%% Optical beerlambert

% Gentot = beerlambert(par, x, source_type, laserlambda, figson)
par_om2.side = 'left';
beerlambert(par_om2, par_om2.xx, 'AM15', 0, true);
beerlambert(par_om2, par_om2.xx, 'laser', 500, true);
par_om2.side = 'right';
beerlambert(par_om2, par_om2.xx, 'AM15', 0, true);
beerlambert(par_om2, par_om2.xx, 'laser', 500, true);

%% Optical lightsource

% AM15 = lightsource(source_type, lambda)
lightsource('AM15', 500);

%% Optical LoadRefrIndex

% [n_interp, k_interp] = LoadRefrIndex(name,wavelengths)
LoadRefrIndex(par.material{1},300:767);

%% Input_files
inputs = dir('Input_files');
for i=1:length(inputs)
    input = inputs(i).name;
    disp(input);
    if ~any(regexp(input,'^\.'))
        par = pc(input);
        xpoints = round(1 + sum(par.layer_points));

        soleq = equilibrate(par);
        
        el = soleq.el;
        el_s = size(el.u);

        exp_el_s = [el.par.tpoints, xpoints, 3];
        assert(all(el_s == exp_el_s), [input ': Expected size: %d, %d, %d. Obtained size: %d, %d, %d.'], el_s(1), el_s(2), el_s(3), exp_el_s(1), exp_el_s(2), exp_el_s(3));
        assert(~any(isnan(el.u(:))))
        
        if par.N_ionic_species > 0
            ion = soleq.ion;
            ion_s = size(ion.u);
            
            exp_ion_s = [ion.par.tpoints, xpoints, round(3+par.N_ionic_species)];
            assert(all(ion_s == exp_ion_s), [input ': Expected size: %d, %d, %d. Obtained size: %d, %d, %d.'], el_s(1), el_s(2), el_s(3), exp_el_s(1), exp_el_s(2), exp_el_s(3));
            assert(~any(isnan(ion.u(:))))
        end
    end
end

%% Scripts
inputs = dir('Scripts');
for i = 1:length(inputs)
    input = inputs(i).name;
    if ~any(regexp(input,'^\.|^test'))
        
        disp(['### running script ' input]);
        run(input);
    end
end

%%
exsol = parex_dactive_light;
%% Helper explore plotPL

% plotPL(exsol)
explore.plotPL(exsol);

%% Helper explore plotsurf

% plotsurf(exsol, yproperty, xlogon, ylogon, zlogon)
explore.plotsurf(exsol, 'Voc_r', true, false, false)

%% Helper explore getJtot

% Jtot = getJtot(sol)
explore.getJtot(soleq.ion);

%% Helper explore writevar

% var = writevar(var, i, j, xx, arr)
x=0;
explore.writevar(x, 1, 1, par.xx, 1);

%% Helper explore helper

% par = helper(par, parname, parvalue)
explore.helper(par, 'Rs', 1);

%% Helper explore plotstat_2D_parval1

% plotstat_2D_parval1(exsol, yproperty, logx, logy)
explore.plotstat_2D_parval1(exsol, 'Voc_r', true, false)

%% Helper explore plotstat_2D_parval2

% plotstat_2D_parval2(exsol, yproperty, logx, logy)
explore.plotstat_2D_parval2(exsol, 'Voc_r', true, false)

%% Helper explore plotfinalELx

% EXPLORE contains functions that plot the final time point solution but
% this functionality is currently not working.
% % plotfinalELx(exsol)
% explore.plotfinalELx(exsol)
% 
% %% Helper explore plotprof_2D
% 
% % plotprof_2D(exsol, yproperty, par1logical, par2logical, logx,logy)
% explore.plotprof_2D(exsol, 'J_f', [true,true], [true,true], true, true)
% 
% %% Helper explore plotU
% 
% % plotU(exsol, par1logical, par2logical,logx,logy)
% explore.plotU(exsol, [true,true], [true,true],false,false)
% 
% %% Helper explore plotCE
% 
% % plotCE(exsol_Voc, exsol_eq, xlogon, ylogon, zlogon, normalise)
% explore.plotCE(exsol, exsol, false, false, false, "ciaomamma")

%% Helper explore plotJV

% plotJV(exsol, par1logical, par2logical)
explore.plotJV(exsol, [true, true, false, true], [true, false, true])

%------------- END OF CODE --------------

