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
par = pc('Input_files/spiro_mapi_tio2.csv');
soleq = equilibrate(par);
solJV = doJV(soleq.ion, 1e-2, 100, 1, 1, 0, 1.4, 3);
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

%% dfplot Jt(sol, xpos)

dfplot.Jt(solJV.ill.f, 2e-5)

%% dfplot Jx(varargin)

dfplot.Jx(soleq.ion)

%% dfplot jx(varargin)

dfplot.jx(soleq.ion)

%% dfplot JV(JV, option)

dfplot.JV(solJV, 3)

%% dfplot Jddx(varargin)

dfplot.Jddx(soleq.ion)

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

%% dfplot npx(varargin)

dfplot.npx(soleq.ion)

%% dfplot acx(varargin)

dfplot.acx(soleq.ion)

%% dfplot gx(varargin)

dfplot.gx(solJV.ill.f)

%% dfplot gxt(sol)

dfplot.gxt(solJV.ill.f)

%% dfplot rx(varargin)

dfplot.rx(solJV.ill.f)

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

%% dfplot deltarhox(varargin)

dfplot.deltarhox(solJV.ill.f)

%% dfplot rhoxFxVx(varargin)

dfplot.rhoxFxVx(soleq.ion)


%% dfplot rhoxVx(varargin)

dfplot.rhoxVx(soleq.ion)

%% dfplot ELx_single(varargin)

dfplot.ELx_single(soleq.ion)

%% dfplot ELnpx(varargin)

dfplot.ELnpx(soleq.ion)

%% dfplot Vacx(varargin)

dfplot.Vacx(soleq.ion)

%% dfplot Vionacx(varargin)

dfplot.Vionacx(soleq.ion)

%% dfplot Fiont(sol, xpos)

dfplot.Fiont(solJV.ill.f, 2e-5)

%% dfplot PLx(varargin)

dfplot.PLx(solJV.ill.f)

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
for i=1:10
   a=[a,0]; 
end

%% Helper calcJsc and calcR0
% [EgArr, Jsc_vs_Eg] = calcJsc
[EgArr, Jsc_vs_Eg] = calcJsc;

% [JV_ana, r0, k_rad, Voc, g0] = calcR0(EgArr, Jsc_vs_Eg, par)
[JV_ana, r0, k_rad, Voc, g0] = calcR0(EgArr, Jsc_vs_Eg, par);

%% Helper Eg_vs_Voc
% [G0_Arr, k_rad_Arr, R0_Arr, Eg, VocArr, DeltaVoc] = Eg_vs_Voc(EgArr, Jsc_vs_Eg)
Eg_vs_Voc

%% Helper explore explore2par
% exsol = explore2par(par_base, parnames, parvalues, JVpnts)
explore.explore2par

%% Helper explore getJtot
% Jtot = getJtot(sol)
getJtot

%% Helper explore writevar
% var = writevar(var, i, j, xx, arr)
writevar

%% Helper explore helper
% par = helper(par, parname, parvalue)
helper

%% Helper explore plotPL
% plotPL(exsol)
plotPL

%% Helper explore plotsurf
% plotsurf(exsol, yproperty, xlogon, ylogon, zlogon)
plotsurf

%% Helper explore plotstat_2D_parval1
% plotstat_2D_parval1(exsol, yproperty, logx, logy)
plotstat_2D_parval1

%% Helper explore plotstat_2D_parval2
% plotstat_2D_parval2(exsol, yproperty, logx, logy)
plotstat_2D_parval2

%% Helper explore plotfinalELx
% plotfinalELx(exsol)
plotfinalELx

%% Helper explore plotprof_2D
% plotprof_2D(exsol, yproperty, par1logical, par2logical, logx,logy)
plotprof_2D

%% Helper explore plotU
% plotU(exsol, par1logical, par2logical,logx,logy)
plotU

%% Helper explore plotCE
% plotCE(exsol_Voc, exsol_eq, xlogon, ylogon, zlogon, normalise)
plotCE

%% Helper explore plotJV
% plotJV(exsol, par1logical, par2logical)
plotJV        
            
%% Helper getpointpos
% ppos = getpointpos(xpos, xmesh)
getpointpos

%% Helper importsolcoresol
% solstruct = importsolcoresol(varargin)
importsolcoresol

%% Helper makemovie
% Framefile = makemovie(sol, plotfun, xrange, yrange, movie_name, Vcounter, tcounter)
makemovie

%% Helper verifyStabilization
% all_stable = verifyStabilization(sol_matrix, t_array, time_fraction)
verifyStabilization

%------------- END OF CODE --------------

