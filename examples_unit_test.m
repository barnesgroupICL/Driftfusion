%EXAMPLES_UNIT_TEST - Runs various simulations in order to test the functioning of pindrift and many other implemented functions
% For profiling the time spent running each function run 'profile on' before running the tests and
% 'profile viewer' after.
% Using Matlab's Coverage Reports, the unused code can be easily spotted.
%
% Syntax:  runtests('examples_unit_test')
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

input_csv = 'Input_files/unit_test_ptpd_mapi_pcbm.csv';

% par = pc(varargin)
par = pc(input_csv);

soleq = equilibrate(par);
soleq_el_1V = jumptoV(soleq.el, 1, 1, false, 0, true, true);
soleq_ion_1V = jumptoV(soleq.ion, 1, 1, true, 0, true, true);
soleq_el_1sun_oc = findVoc(soleq.el, 1, false, 0.9, 1.3, 1e-5);
soleq_ion_1sun_oc = findVoc(soleq.ion, 1, true, 0.9, 1.3, 1e-5);
structs_el_sc_ill = genIntStructs(soleq.el, 1, 1e-3, 2, true);
structs_ion_sc_ill = genIntStructs(soleq.ion, 1, 1e-3, 2, true);
structs_el_vapp = genVappStructs(soleq.el, [0.6,0.9,1.2], false);

el_dark_0V_1Hz = doIS_EA(soleq.el, 1e-3, 1, 5, 8, 1, 1e-2);
ion_dark_0V_1Hz = doIS_EA(soleq.ion, 1e-3, 1, 5, 8, 1, 1e-2);

structs_ion_sc_ill_IS = IS_script(structs_ion_sc_ill, 1e6, 1e-3, 4, 1e-3, false, true, false);
structs_el_vapp_IS = IS_script(structs_el_vapp, 1e6, 1e-3, 4, 1e-3, false, true, false);

structs_ion_sc_ill_EA = EA_script(structs_ion_sc_ill, 1e6, 1e-3, 4, 1e-3, false, false);
structs_el_vapp_EA = EA_script(structs_el_vapp, 1e6, 1e-3, 4, 1e-3, false, false);

t = (0:0.1:12)';
sin_func = fun_gen('sin');
coeff1 = [123456, 0.123456, 1, deg2rad(0.01)];
y1 = sin_func(coeff1, t);
coeff2 = [-123456, 0.123456, 1, deg2rad(-0.01)];
y2 = sin_func(coeff2, t);
coeff3 = [0.123456, 123456, 1, deg2rad(89.999)];
y3 = sin_func(coeff3, t);
coeff4 = [-0.123456, 123456, 1, deg2rad(-89.999)];
y4 = sin_func(coeff4, t);


%% Protocols doIS_EA el dark 0V lowfreq reachStability
% struct_IS = doIS_EA(struct_Int, deltaV, freq, periods, tpoints_per_period, stability_timefraction, RelTol)

doIS_EA(soleq.el, 1e-3, 1e-2, 5, 8, 0.5, 1e-4);

%% Protocols doIS_EA el dark 0V lowfreq smallOscillations
doIS_EA(soleq.el, 1e-3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA el dark 0V lowfreq largeOscillations
doIS_EA(soleq.el, 0.3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA el dark 0V highfreq
doIS_EA(soleq.el, 1e-3, 1e7, 5, 8, 1, 1e-2);



%% Protocols doIS_EA el dark 1V lowfreq reachStability
doIS_EA(soleq_el_1V, 1e-3, 1e-2, 5, 8, 0.5, 1e-4);

%% Protocols doIS_EA el dark 1V lowfreq smallOscillations
doIS_EA(soleq_el_1V, 1e-3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA el dark 1V lowfreq largeOscillations
doIS_EA(soleq_el_1V, 0.3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA el dark 1V highfreq
doIS_EA(soleq_el_1V, 1e-3, 1e7, 5, 8, 1, 1e-2);



%% Protocols doIS_EA el 1sun OC lowfreq reachStability
doIS_EA(soleq_el_1sun_oc, 1e-3, 1e-2, 5, 8, 0.5, 1e-4);

%% Protocols doIS_EA el 1sun OC lowfreq smallOscillations
doIS_EA(soleq_el_1sun_oc, 1e-3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA el 1sun OC lowfreq largeOscillations
doIS_EA(soleq_el_1sun_oc, 0.3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA el 1sun OC highfreq
doIS_EA(soleq_el_1sun_oc, 1e-3, 1e7, 5, 8, 1, 1e-2);



%% Protocols doIS_EA el 1sun OC lowfreq reachStability
doIS_EA(structs_el_sc_ill{1,end}, 1e-3, 1e-2, 5, 8, 0.5, 1e-4);

%% Protocols doIS_EA el 1sun OC lowfreq smallOscillations
doIS_EA(structs_el_sc_ill{1,end}, 1e-3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA el 1sun OC lowfreq largeOscillations
doIS_EA(structs_el_sc_ill{1,end}, 0.3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA el 1sun OC highfreq
doIS_EA(structs_el_sc_ill{1,end}, 1e-3, 1e7, 5, 8, 1, 1e-2);



%% Protocols doIS_EA ion dark 0V lowfreq reachStability
doIS_EA(soleq.ion, 1e-3, 1e-2, 5, 8, 0.5, 1e-4);

%% Protocols doIS_EA ion dark 0V lowfreq smallOscillations
doIS_EA(soleq.ion, 1e-3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA ion dark 0V lowfreq largeOscillations
doIS_EA(soleq.ion, 0.3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA ion dark 0V highfreq
doIS_EA(soleq.ion, 1e-3, 1e7, 5, 8, 1, 1e-2);



%% Protocols doIS_EA ion dark 1V lowfreq reachStability
doIS_EA(soleq_ion_1V, 1e-3, 1e-2, 5, 8, 0.5, 1e-4);

%% Protocols doIS_EA ion dark 1V lowfreq smallOscillations
doIS_EA(soleq_ion_1V, 1e-3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA ion dark 1V lowfreq largeOscillations
doIS_EA(soleq_ion_1V, 0.3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA ion dark 1V highfreq
doIS_EA(soleq_ion_1V, 1e-3, 1e7, 5, 8, 1, 1e-2);



%% Protocols doIS_EA ion 1sun OC lowfreq reachStability
doIS_EA(soleq_ion_1sun_oc, 1e-3, 1e-2, 5, 8, 0.5, 1e-4);

%% Protocols doIS_EA ion 1sun OC lowfreq smallOscillations
doIS_EA(soleq_ion_1sun_oc, 1e-3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA ion 1sun OC lowfreq largeOscillations
doIS_EA(soleq_ion_1sun_oc, 0.3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA ion 1sun OC highfreq
doIS_EA(soleq_ion_1sun_oc, 1e-3, 1e7, 5, 8, 1, 1e-2);



%% Protocols doIS_EA ion 1sun SC lowfreq reachStability
doIS_EA(structs_ion_sc_ill{1,end}, 1e-3, 1e-2, 5, 8, 0.5, 1e-4);

%% Protocols doIS_EA ion 1sun SC lowfreq smallOscillations
doIS_EA(structs_ion_sc_ill{1,end}, 1e-3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA ion 1sun SC lowfreq largeOscillations
doIS_EA(structs_ion_sc_ill{1,end}, 0.3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA ion 1sun SC highfreq
doIS_EA(structs_ion_sc_ill{1,end}, 1e-3, 1e7, 5, 8, 1, 1e-2);



%% Protocols doIS_EA el dark Vapp lowfreq reachStability
doIS_EA(structs_el_vapp{1,end}, 1e-3, 1e-2, 5, 8, 0.5, 1e-4);

%% Protocols doIS_EA el dark Vapp lowfreq smallOscillations
doIS_EA(structs_el_vapp{1,end}, 1e-3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA el dark Vapp lowfreq largeOscillations
doIS_EA(structs_el_vapp{1,end}, 0.3, 1e-2, 5, 8, 1, 1e-2);

%% Protocols doIS_EA el dark Vapp highfreq
doIS_EA(structs_el_vapp{1,end}, 1e-3, 1e7, 5, 8, 1, 1e-2);



%% Analysis IS_EA_ana_demodulation
% coeff = IS_EA_ana_demodulation(t, y, fun_type, freq)

coeff = IS_EA_ana_demodulation(t, y1, 'sin', 1);
assert(ismembertol(coeff, coeff1([1,2,4]), 1e-6, 'ByRows', true))

coeff = IS_EA_ana_demodulation(t, y2, 'sin', 1);
assert(ismembertol(coeff, coeff2([1,2,4]), 1e-6, 'ByRows', true))

coeff = IS_EA_ana_demodulation(t, y3, 'sin', 1);
assert(ismembertol(coeff, coeff3([1,2,4]), 1e-6, 'ByRows', true))

coeff = IS_EA_ana_demodulation(t, y4, 'sin', 1);
assert(ismembertol(coeff, coeff4([1,2,4]), 1e-6, 'ByRows', true))



%% Analysis IS_EA_ana_fit
% coeff = IS_EA_ana_fit(t, y, fun_type, freq)

coeff = IS_EA_ana_fit(t, y1, 'sin', 1);
assert(ismembertol(coeff, coeff1([1,2,4]), 1e-6, 'ByRows', true))

coeff = IS_EA_ana_fit(t, y2, 'sin', 1);
assert(ismembertol(coeff, coeff2([1,2,4]), 1e-6, 'ByRows', true))

coeff = IS_EA_ana_fit(t, y3, 'sin', 1);
assert(ismembertol(coeff, coeff3([1,2,4]), 1e-6, 'ByRows', true))

coeff = IS_EA_ana_fit(t, y4, 'sin', 1);
assert(ismembertol(coeff, coeff4([1,2,4]), 1e-6, 'ByRows', true))



%% Analysis IS_ana
% [n_coeff, i_coeff, r_coeff, dQ_coeff, n_noionic_coeff] = IS_ana(struct_IS, minimal_mode, demodulation)

IS_ana(el_dark_0V_1Hz, false, true);

IS_ana(ion_dark_0V_1Hz, false, false);

IS_ana(ion_dark_0V_1Hz, true, true);



%% Analysis IS_ana_subtracting
% [subtracting_q_t, subtracting_q_intr_t] = IS_ana_subtracting(struct_IS)

IS_ana_subtracting(ion_dark_0V_1Hz);



%% Analysis EA_ana
% EA_ana(struct_IS_EA, do_graphics, local_field, demodulation, savefig_dir)

EA_ana(el_dark_0V_1Hz, true, true, true, missing);

EA_ana(ion_dark_0V_1Hz, false, true, false, missing);

EA_ana(ion_dark_0V_1Hz, true, false, false, 'test_EAana_delete_me');



%% Scripts IS_script_nonparallel el SC ill
% IS_results = IS_script_nonparallel(structs, startFreq, endFreq, Freq_points, deltaV, sequential, frozen_ions, demodulation, do_graphics, save_solutions)

IS_script_nonparallel(structs_el_sc_ill, 1e9, 1e-3, 3, 1e-3, false, false, true, false, false);

%% Scripts IS_script_nonparallel el Vapp
IS_script_nonparallel(structs_el_vapp, 1e9, 1e-3, 3, 1e-3, false, false, true, false, false);

%% Scripts IS_script_nonparallel ion SC ill
IS_script_nonparallel(structs_ion_sc_ill, 1e9, 1e-3, 3, 1e-3, false, false, true, false, false);

%% Scripts IS_script_nonparallel el SC dark
IS_script_nonparallel(soleq.el, 1e9, 1e-3, 3, 1e-3, false, false, true, false, false);

%% Scripts IS_script_nonparallel ion SC dark sequential
IS_script_nonparallel(soleq.ion, 1e9, 1e-1, 6, 1e-3, true, false, true, false, false);

%% Scripts IS_script_nonparallel ion SC dark frozenIons
IS_script_nonparallel(soleq.ion, 1e3, 1e-3, 3, 1e-3, false, true, true, false, false);

%% Scripts IS_script_nonparallel el SC dark fit
IS_script_nonparallel(soleq.el, 1e9, 1e-3, 3, 1e-3, false, false, false, false, false);

%% Scripts IS_script_nonparallel el SC dark graphics
IS_script_nonparallel(soleq.el, 1e9, 1e-3, 3, 1e-3, false, false, true, true, false);

%% Scripts IS_script_nonparallel el SC dark savesolutions
IS_script_nonparallel(soleq.el, 1e9, 1e-3, 3, 1e-3, false, false, true, false, true);



%% Scripts IS_script el SC ill
% IS_results = IS_script(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, demodulation, do_graphics)

IS_script(structs_el_sc_ill, 1e9, 1e-3, 3, 1e-3, false, true, false);

%% Scripts IS_script el Vapp
IS_script(structs_el_vapp, 1e9, 1e-3, 3, 1e-3, false, true, false);

%% Scripts IS_script ion SC ill
IS_script(structs_ion_sc_ill, 1e9, 1e-3, 3, 1e-3, false, true, false);

%% Scripts IS_script el SC dark
IS_script(soleq.el, 1e9, 1e-3, 3, 1e-3, false, true, false);

%% Scripts IS_script ion SC dark 
IS_script(soleq.ion, 1e9, 1e-1, 6, 1e-3, true, true, false);

%% Scripts IS_script ion SC dark frozenIons
IS_script(soleq.ion, 1e3, 1e-3, 3, 1e-3, true, true, false);

%% Scripts IS_script el SC dark fit
IS_script(soleq.el, 1e9, 1e-3, 3, 1e-3, false, false, false);

%% Scripts IS_script el SC dark graphics
IS_script(soleq.el, 1e9, 1e-3, 3, 1e-3, false, true, true);



%% Scripts EA_script el SC ill
% EA_results = EA_script(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, do_graphics)

EA_script(structs_el_sc_ill, 1e9, 1e-3, 3, 1e-3, false, false);

%% Scripts EA_script el Vapp
EA_script(structs_el_vapp, 1e9, 1e-3, 3, 1e-3, false, false);

%% Scripts EA_script ion SC ill
EA_script(structs_ion_sc_ill, 1e9, 1e-3, 3, 1e-3, false, false);

%% Scripts EA_script el SC dark
EA_script(soleq.el, 1e9, 1e-3, 3, 1e-3, false, false);

%% Scripts EA_script ion SC dark frozenIons
EA_script(soleq.ion, 1e3, 1e-3, 3, 1e-3, true, false);

%% Scripts EA_script el SC dark graphics
EA_script(soleq.el, 1e9, 1e-3, 3, 1e-3, false, true);



%% Analysis IS_script_ana_impedance
% IS_script_ana_impedance(IS_results)

IS_script_ana_impedance(structs_el_vapp_IS)

IS_script_ana_impedance(structs_ion_sc_ill_IS)



%% Analysis IS_script_ana_nyquist
% IS_script_ana_nyquist(IS_results)

IS_script_ana_nyquist(structs_el_vapp_IS)

IS_script_ana_nyquist(structs_ion_sc_ill_IS)



%% Analysis IS_script_ana_phase
% IS_script_ana_phase(IS_results)

IS_script_ana_phase(structs_el_vapp_IS)

IS_script_ana_phase(structs_ion_sc_ill_IS)



%% Analysis EA_script_ana_Efield
% EA_script_ana_Efield(EA_results)

EA_script_ana_Efield(structs_el_vapp_EA)

EA_script_ana_Efield(structs_ion_sc_ill_EA)



%% Analysis EA_script_ana_phase
% EA_script_ana_phase(EA_results)

EA_script_ana_phase(structs_el_vapp_EA)

EA_script_ana_phase(structs_ion_sc_ill_EA)



%% Helper export_IS_EA_struct
% export_IS_EA_struct(struct, prefix)

export_IS_EA_struct(el_dark_0V_1Hz, 'test_exportISEAstruct_delete_me')

export_IS_EA_struct(ion_dark_0V_1Hz, 'test_exportISEAstruct_delete_me')



%% Helper export_IS_results
% export_IS_results(IS_results, prefix)

export_IS_results(structs_el_vapp_IS, 'test_exportISresults_vapp_delete_me')

export_IS_results(structs_ion_sc_ill_IS, 'test_exportISresults_ill_delete_me')



%------------- END OF CODE --------------

