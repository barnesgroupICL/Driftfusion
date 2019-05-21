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
[sol_eq, sol_i_eq, sol_i_eq_SR, ssol_i_eq, ssol_i_eq_SR, sol_i_1S_SR, ssol_i_1S_SR] = equilibrate_minimal();
structs_oc = genIntStructs(ssol_i_eq_SR, 1, 1e-3, 2, true);
structs_sc = genIntStructs(sol_i_eq_SR, 1, 1e-3, 2, true);
structs_vapp = genVappStructs(sol_i_eq_SR, [0.4, 0.95]);
structs_oc_noions_noSR = genIntStructs(symmetricize(sol_eq), 1, 1e-3, 2, true);
structs_sc_noions_noSR = genIntStructs(sol_eq, 1, 1e-3, 2, true);
[~, ~, sol_i_eq_SR_10kxSRH_001xmajority, ~, ~, ~, ~] = equilibrate_minimal(pinParams_10kxSRH_001xmajority);
structs_vapp_10kxSRH_001xmajority = genVappStructs(sol_i_eq_SR_10kxSRH_001xmajority, [0.4, 0.7]);

% test Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel)
% ISwave_results =  ISwave_full_exec_nonparallel(structs, startFreq, endFreq, Freq_points, deltaV, sequential, frozen_ions, demodulation, do_graphics, save_solutions)

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - OC and various light intensities
results = ISwave_full_exec_nonparallel(structs_oc, 1e2, 1e2, 1, 0.002, false, false, true, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, 1S and various frequencies
results = ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 6, 0.002, false, false, true, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, dark and various frequencies
results = ISwave_full_exec_nonparallel(ssol_i_eq_SR, 1e8, 1e-2, 6, 0.002, false, false, true, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - SC and various light intensities
results = ISwave_full_exec_nonparallel(structs_sc, 1e2, 1e2, 1, 0.002, false, false, true, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, 1S and big oscillations
results = ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e2, 1e2, 1, 0.2, false, false, true, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, 1S and do not reach stability and do sequential simulations
results = ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e6, 1e-2, 5, 0.002, true, false, true, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - do graphics
results = ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, false, false, true, true, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - frozen ions & do graphics & saving solutions
results = ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, false, true, true, true, true);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - applied voltage
results = ISwave_full_exec_nonparallel(structs_vapp, 1e2, 1e2, 1, 0.002, false, false, true, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, no ions and no surface recombination
results = ISwave_full_exec_nonparallel(structs_oc_noions_noSR, 1e2, 1e2, 1, 0.002, false, false, true, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - SC, no ions and no surface recombination
results = ISwave_full_exec_nonparallel(structs_sc_noions_noSR, 1e2, 1e2, 1, 0.002, false, false, true, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - fitting instead of demodulation
results = ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, false, false, false, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage and without parallelization (ISwave_nonparallel) - negative capacitance
results = ISwave_full_exec_nonparallel(structs_vapp_10kxSRH_001xmajority, 1e-2, 1e-2, 1, 0.002, false, false, true, false, false);
assert(~any(isnan(results.J_phase(:))))

% test Impedance Spectroscopy with oscillating voltage (ISwave)
% ISwave_results = ISwave_full_exec(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, demodulation, do_graphics)

%% IS with oscillating voltage (ISwave) - OC and various light intensities
results = ISwave_full_exec(structs_oc, 1e2, 1e2, 1, 0.002, false, true, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage (ISwave) - OC, dark and various frequencies
results = ISwave_full_exec(ssol_i_eq_SR, 1e8, 1e-2, 6, 0.002, false, true, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage (ISwave) - SC, 1S and one frequency
results = ISwave_full_exec(sol_i_1S_SR, 1e2, 1e2, 1, 0.002, false, true, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage (ISwave) - Vapp, do graphics
results = ISwave_full_exec(structs_vapp, 1e8, 1e-2, 2, 0.002, true, true, true);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage (ISwave) - OC, frozen ions & do graphics
results = ISwave_full_exec(structs_oc, 1e8, 1e-2, 2, 0.002, true, true, true);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage (ISwave) - applied voltage
results = ISwave_full_exec(structs_vapp, 1e2, 1e2, 1, 0.002, false, true, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage (ISwave) - no ions and no surface recombination
results = ISwave_full_exec(structs_sc_noions_noSR, 1e2, 1e2, 1, 0.002, false, true, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage (ISwave) - fitting instead of demodulation
results = ISwave_full_exec(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, false, false, false);
assert(~any(isnan(results.J_phase(:))))

%% IS with oscillating voltage (ISwave) - negative capacitance
results = ISwave_full_exec(structs_vapp_10kxSRH_001xmajority, 1e-2, 1e-2, 1, 0.002, false, true, false);
assert(~any(isnan(results.J_phase(:))))

% test Ideality Factor from VOC
% nid = ideality_from_Voc(structs_oc)

%% Ideality Factor from VOC
ideality_from_Voc(structs_oc)

% test Ideality Factor from stabilized JV points in dark
% nid_array = ideality_from_dark_jVpoints(asymstruct_eq, Vend, deltaV)

%% Ideality Factor from stabilized JV points in dark
ideality_from_dark_jVpoints(sol_i_eq_SR, 1, 0.1)

%------------- END OF CODE --------------

