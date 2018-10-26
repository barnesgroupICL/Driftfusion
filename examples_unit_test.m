% test pindrift and photophysics examples
% USAGE: runtests('examples_unit_test')
% for profiling the usage run 'profile on' before running the tests and
% 'profile viewer' after

% prepare solutions and indirectly test equilibrate and genIntStructs
% functions
[sol_eq, sol_i_eq, sol_i_eq_SR, ssol_i_eq, ssol_i_eq_SR, sol_i_1S_SR, ssol_i_1S_SR] = equilibrate_minimal();
symStructs = genIntStructs(ssol_i_eq_SR, 10, 1e-5, 4, true);
asymStructs = genIntStructs(sol_i_eq_SR, 1, 1e-6, 4, true);
vappStructs = genVappStructs(sol_i_eq_SR, [0.2, 0.4, 0.95]);
asymStructsNoIonsNoSR = genIntStructs(sol_eq, 10, 1e-5, 4, true);

% test Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel)
% ISwave_results =  ISwave_full_exec_nonparallel(structs, startFreq, endFreq, Freq_points, deltaV, sequential, frozen_ions, do_graphics, save_solutions)

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC and various light intensities
ISwave_full_exec_nonparallel(symStructs, 1e2, 1e2, 1, 0.002, false, false, true, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, 1S and various frequencies
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 6, 0.002, false, false, true, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, dark and various frequencies
ISwave_full_exec_nonparallel(ssol_i_eq_SR, 1e8, 1e-2, 6, 0.002, false, false, true, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - SC and various light intensities
ISwave_full_exec_nonparallel(asymStructs, 1e2, 1e2, 1, 0.002, false, false, true, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, 1S and big oscillations
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e2, 1e2, 1, 0.2, false, false, true, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, 1S and do not reach stability and do sequential simulations
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 3, 0.002, true, false, true, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - do graphics
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, false, false, true, true, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - frozen ions & do graphics & saving solutions
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, false, true, true, true, true);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - applied voltage
ISwave_full_exec_nonparallel(vappStructs, 1e2, 1e2, 1, 0.002, false, false, true, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - no ions and no surface recombination
ISwave_full_exec_nonparallel(asymStructsNoIonsNoSR, 1e2, 1e2, 1, 0.002, false, false, true, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - fitting instead of demodulation
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, false, false, false, false, false);

% test Impedance Spectroscopy with oscillating voltage (ISwave)
% ISwave_results = ISwave_full_exec(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, do_graphics)

%% Impedance Spectroscopy with oscillating voltage (ISwave) - OC and various light intensities
ISwave_full_exec(symStructs, 1e2, 1e2, 1, 0.002, false, true, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - OC, dark and various frequencies
ISwave_full_exec(ssol_i_eq_SR, 1e8, 1e-2, 6, 0.002, false, true, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - SC, 1S and one frequency
ISwave_full_exec(sol_i_1S_SR, 1e2, 1e2, 1, 0.002, false, true, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - frozen ions & do graphics
ISwave_full_exec(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, true, true, true);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - applied voltage
ISwave_full_exec(vappStructs, 1e2, 1e2, 1, 0.002, false, true, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - no ions and no surface recombination
ISwave_full_exec(asymStructsNoIonsNoSR, 1e2, 1e2, 1, 0.002, false, true, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - fitting instead of demodulation
ISwave_full_exec(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, false, false, false);

