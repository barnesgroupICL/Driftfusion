%EXAMPLES_UNIT_TEST - Runs various simulations in order to test the functioning of pindrift and many other implemented functions
% For profiling the time spent running each function run 'profile on' before running the tests and
% 'profile viewer' after.
% Using Matlab's Coverage Reports, the unused code can be easily spotted.
%
% Syntax:  runtests('unit_test')
%
% Other m-files required: pindrift
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: May 2018

%------------- BEGIN CODE --------------

% prepare solutions and indirectly test equilibrate and genIntStructs
% functions
initialise_df
par = pc('Input_files/ptpd_mapi_pcbm.csv');
soleq = equilibrate(par);
solJV = doJV(soleq.ion, 1e-2, 50, 1, 1, 0, 1.2, 3);
% example taken from Scripts/explore_script
exsol = explore.explore2par(par, {'d(1,3)','Int'},...
    {[400e-7, 800e-7], logspace(-2,0,2)}, 20);
% structs_oc = genIntStructs(ssol_i_eq_SR, 1, 1e-3, 2, true);
% structs_sc = genIntStructs(sol_i_eq_SR, 1, 1e-3, 2, true);
% structs_vapp = genVappStructs(sol_i_eq_SR, [0.4, 0.95]);
% structs_oc_noions_noSR = genIntStructs(symmetricize(sol_eq), 1, 1e-3, 2, true);
% structs_sc_noions_noSR = genIntStructs(sol_eq, 1, 1e-3, 2, true);
% [~, ~, sol_i_eq_SR_10kxSRH_001xmajority, ~, ~, ~, ~] = equilibrate_minimal(pinParams_10kxSRH_001xmajority);
% structs_vapp_10kxSRH_001xmajority = genVappStructs(sol_i_eq_SR_10kxSRH_001xmajority, [0.4, 0.7]);

% test Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel)
% ISwave_results =  ISwave_full_exec_nonparallel(structs, startFreq, endFreq, Freq_points, deltaV, sequential, frozen_ions, demodulation, do_graphics, save_solutions)

%% df no input

df();
 
%% df one input

df(soleq.ion);

%% df two inputs zero u

sol = struct();
sol.u = 0;
df(sol, par);

%% df two inputs char

df(soleq.ion, 'test');

%% df two inputs else

df(soleq.ion, par);

%% creation equilibrium solutions

inputs = dir('Input_files');
for i=1:length(inputs)
    input = inputs(i).name;
    if ~any(regexp(input,'^\.'))
        
        disp(['### ' input ' parameters generation']);
        par = pc(input);
        xpoints = round(1 + sum(par.layer_points));

        disp(['### ' input ' generation of solutions at equilibrium']);
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

%% example scripts

inputs = dir('Scripts');
for i=1:length(inputs)
    input = inputs(i).name;
    if ~any(regexp(input,'^\.|^test'))
        
        disp(['### running script ' input]);
        run(input);
    end
end

%% dfplot ELx(varargin)

dfplot.ELx(soleq.ion)
dfplot.ELx(JVsol.ill.f,[0,100])


%% dfplot Jt(sol, xpos)

dfplot.Jt(solJV.ill.f, 2e-5)

%% dfplot Jx(varargin)

dfplot.Jx(soleq.ion)
dfplot.Jx(JVsol.ill.f,[0,100])

%% dfplot jx(varargin)

dfplot.jx(soleq.ion)
dfplot.jx(JVsol.ill.f,[0,100])

%% dfplot JV(JV, option)

dfplot.JV(solJV, 3)

%% dfplot Jddx(varargin)

dfplot.Jddx(soleq.ion)
dfplot.Jddx(JVsol.ill.f,[0,100])

%% dfplot Voct(sol)

dfplot.Voct(solJV.ill.f)

%% dfplot PLt(sol)

dfplot.PLt(solJV.ill.f)

%% dfplot Vappt(sol)

dfplot.Vappt(solJV.ill.f)

%% dfplot JVapp(sol, xpos)

dfplot.JVapp(soleq.ion, 2e-5)

%% dfplot JtotVapp(sol, xpos)

dfplot.JtotVapp(solJV.ill.f, 2e-5)

%% dfplot logJVapp(sol, xpos)

dfplot.logJVapp(solJV.ill.f, 2e-5)

%% dfplot logJVapp3D(sol, xpos, ylogon)

dfplot.logJVapp3D(solJV.ill.f, 2e-5, true)

%% dfplot xmesh(sol)

dfplot.xmesh(soleq.ion)

%% dfplot Vx(varargin)

dfplot.Vx(soleq.ion)
dfplot.Vx(JVsol.ill.f,[0,100])

%% dfplot npx(varargin)

dfplot.npx(soleq.ion)
dfplot.npx(JVsol.ill.f,[0,100])

%% dfplot acx(varargin)

dfplot.acx(soleq.ion)
dfplot.acx(JVsol.ill.f,[0,100])

%% dfplot gx(varargin)

dfplot.gx(solJV.ill.f)
dfplot.gx(JVsol.ill.f,[0,100])

%% dfplot gxt(sol)

dfplot.gxt(solJV.ill.f)

%% dfplot rx(varargin)

dfplot.rx(solJV.ill.f)
dfplot.rx(JVsol.ill.f,[0,100])

%% dfplot JVrec(JV, option)

dfplot.JVrec(solJV, 3)

%% dfplot Ft(sol, xpos)

dfplot.Ft(solJV.ill.f, 2e-5)

%% dfplot sigmat(sol)

dfplot.sigmat(solJV.ill.f)

%% dfplot Qt(sol, x1, x2)

dfplot.Qt(solJV.ill.f)

%% dfplot QVapp(sol, x1, x2)

dfplot.QVapp(solJV.ill.f)

%% dfplot rhox(varargin)

dfplot.rhox(soleq.ion)
dfplot.rhox(JVsol.ill.f,[0,100])

%% dfplot deltarhox(varargin)

dfplot.deltarhox(solJV.ill.f)
dfplot.deltarhox(JVsol.ill.f,[0,100])

%% dfplot rhoxFxVx(varargin)

dfplot.rhoxFxVx(soleq.ion)
dfplot.rhoxFxVx(JVsol.ill.f,[0,100])

%% dfplot rhoxVx(varargin)

dfplot.rhoxVx(soleq.ion)
dfplot.rhoxVx(JVsol.ill.f,[0,100])

%% dfplot ELx_single(varargin)

dfplot.ELx_single(soleq.ion)
dfplot.ELx_single(JVsol.ill.f,[0,100])

%% dfplot ELnpx(varargin)

dfplot.ELnpx(soleq.ion)
dfplot.ELnpx(JVsol.ill.f,[0,100])

%% dfplot Vacx(varargin)

dfplot.Vacx(soleq.ion)
dfplot.Vacx(JVsol.ill.f,[0,100])

%% dfplot Vionacx(varargin)

dfplot.Vionacx(soleq.ion)
dfplot.Vionacx(JVsol.ill.f,[0,100])

%% dfplot Fiont(sol, xpos)

dfplot.Fiont(solJV.ill.f, 2e-5)

%% dfplot PLx(varargin)

dfplot.PLx(solJV.ill.f)
dfplot.PLx(JVsol.ill.f,[0,100])

%% dfplot colourblocks(sol, yrange)

dfplot.colourblocks(soleq.ion, [0,1])

%% dfplot sortarg(args)

test_tarr = [0,10,100];
test_xrange = [150, 250];
[~, ~, ~, ~] = dfplot.sortarg({solJV.ill.f});
[~, tarr, ~, ~] = dfplot.sortarg({solJV.ill.f, test_tarr});
assert(all(tarr == test_tarr))
[~, tarr, ~, xrange] = dfplot.sortarg({solJV.ill.f, test_tarr, test_xrange});
assert(all([tarr == test_tarr, xrange == test_xrange]))

%% dfplot x2d(sol, xmesh, variables, legstr, linestyle, ylab, tarr, xrange, logx, logy)

dfplot.x2d(solJV.ill.f, solJV.ill.f.x, {solJV.ill.f.u(:,:,1)}, {'test'}, ['-','.'], 'ylab', solJV.ill.f.t, [150,250], false, false)

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
Framefile = makemovie(solJV.ill.f, @dfplot.PLx, [0,100e-7], 0, 'test_makemovie_delete_me', true, true);

%% Helper verifyStabilization
% all_stable = verifyStabilization(sol_matrix, t_array, time_fraction)
all_stable = verifyStabilization(soleq.ion.u, soleq.ion.t, 0.1);

%% Protocols changeLight.m
% sol_int = changeLight(sol, newInt, tmax)
sol_int1 = changeLight(soleq.ion, 1, 2);

%% Protocols doCV.m
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV = doCV(soleq.ion, 0.1, 0, 0.5, -0.5, 0.1, 1, 50);

%% Protocols doIMPS.m
% sol_IMPS = doIMPS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
sol_IMPS = doIMPS(soleq.ion, 0.2, 0.2, 1, 10, 50);

%% Protocols doIMVS.m
% sol_IMVS = doIMVS(sol_ini, int_base, int_delta, frequency, tmax, tpoints)
sol_IMVS = doIMVS(soleq.ion, 0.2, 0.2, 1, 10, 50);

%% Protocols doLightPulse.m
% [sol_pulse] = doLightPulse(sol_ini, pulse_int, tmax, tpoints, duty, mobseti, log_timemesh)
sol_pulse = doLightPulse(soleq.ion, 0.1, 1, 20, 5, true, true);
sol_pulse2 = doLightPulse(soleq.ion, 0.1, 1, 20, 5, true, false);

%% Protocols doSDP.m and Analysis anasdp
% sdpsol = doSDP(sol_ini, tdwell_arr, Vjump, pulse_int, pulse_tmax, duty, tpoints, scalefactor)
sdpsol = doSDP(soleq.ion, [1e-3,1], 0.9, 1, 1, 5, 50, 1);

% anasdp(sdpsol, Jtr_time)
anasdp(sdpsol, 1);

%% Protocols doSPV.m and Analysis spvana
% spvsol = doSPV(sol_ini, Int, mobseti, tpoints, tmax, Rs, stabilise)
spvsol = doSPV(soleq.ion, 1, true, 20, 1, 1e3, true);

% spvdat = spvana(spvsol)
spvdat = spvana(spvsol);

%% Protocols doTPV.m
% [sol_TPV, sol_ill] = doTPV(sol_ini, bias_int, stab_time, mobseti, Rs, pulse_int, tmax, tpoints, duty)
[sol_TPV, sol_ill] = doTPV(soleq.ion, 1, 10, true, 1e-2, 0.1, 1, 20, 5);

%% Protocols findVocDirect.m
% [sol_Voc, Voc] = findVocDirect(sol_ini, light_intensity, mobseti)
[sol_Voc, Voc] = findVocDirect(soleq.ion, 1, true);

%% Protocols findVoc.m
% [sol_Voc, Voc] = findVoc(sol_ini, Int, mobseti, x0, x1, tol)
[sol_Voc, Voc] = findVoc(soleq.ion, 1, true, 0.9, 1.3, 1e-5)

%% Protocols genIntStructs.m
% [structCell, V_array, J_array] = genIntStructs(struct_eq, startInt, endInt, points, include_dark)
[structCell, V_array, J_array] = genIntStructs(soleq.ion, 1e-4, 1, 5, true);

%% Protocols genIntStructsRealVoc.m
% [goodVocAsymStructCell, VOCs] = genIntStructsRealVoc(struct_eq, startInt, endInt, points, include_dark)
[goodVocAsymStructCell, VOCs] = genIntStructsRealVoc(soleq.ion, 1e-4, 1e-3, 2, false);

%% Protocols genVappStructs.m
% VappSol = genVappStructs(solini, Vapp_arr, mobseti)
VappSol = genVappStructs(soleq.ion, [0.9,1.1,1.3], true);

%% Protocols jumptoV.m
% sol_relax = jumptoV(sol_ini, Vjump, tdwell, mobseti, Int, stabilise, accelerate)
sol_relax = jumptoV(soleq.ion, 1.2, 1, true, 1, true, true);

%% Protocols lightonRs.m
% sol_ill = lightonRs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
sol_ill = lightonRs(soleq.ion, 1, 10, true, 1e6, 50);

%% Protocols stabilize.m
% steadystate_struct = stabilize(struct)
steadystate_struct = stabilize(struct);

%% Protocols sweeplight.m
% sol_sweep = sweeplight(sol_ini, tmax, tpoints, end_int)
sol_sweep = sweeplight(soleq.ion, 1, 50, 1);

%% Protocols transient_nid.m
% sol_OC = transient_nid(sol_ini, int_arr, stab_time, mobseti, Rs, pnts)

%% Protocols VappFunction.m
% sol = VappFunction(sol_ini, Vapp_func, Vapp_coeff, tmax, tpoints, logtime)

%------------- END OF CODE --------------

