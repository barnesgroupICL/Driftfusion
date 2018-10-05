% test pindrift and photophysics examples
% USAGE: runtests('examples_unit_test')
% for profiling the usage run 'profile on' before running the tests and
% 'profile viewer' after

% prepare solutions and indirectly test equilibrate and genIntStructs
% functions
[sol_eq, sol_i_eq, sol_i_eq_SR, ssol_i_eq, ssol_i_eq_SR, sol_i_1S_SR, ssol_i_1S_SR] = equilibrate_minimal();
symstructs = genIntStructs(ssol_i_eq_SR, 10, 1e-5, 4, true);
asymstructs = genIntStructs(sol_i_eq_SR, 10, 1e-5, 4, true);

% test Charge Extraction (CE)
% CE_struct = CE_full_exec(symstructs, BC, save_solutions, save_results)

%% Charge Extraction (CE) - various light intensities
%CE_full_exec(symstructs, 1, false, false);

%% Charge Extraction (CE) - just dark & saving solutions
%CE_full_exec(ssol_i_eq, 1, true, true);

%% Charge Extraction (CE) - boundary conditions and just one sun
%CE_full_exec(ssol_i_light_BC2, 2, false, false);

% test Transient PhotoVoltage (TPV) with constant pulse intensity
% TPV_struct = TPVconst_full_exec(symstructs, save_solutions, save_results)

%% Transient PhotoVoltage (TPV) - constant pulse intensity, various light intensities
%TPVconst_full_exec(symstructs, false, false);

%% Transient PhotoVoltage (TPV) - constant pulse intensity, just dark & saving solutions
%TPVconst_full_exec(ssol_i_eq, true, true);

%% Transient PhotoVoltage (TPV) - constant pulse intensity, just one sun
%TPVconst_full_exec(ssol_i_light, false, false);

% test Transient PhotoVoltage (TPV) with variable pulse intensity
% TPV_struct = TPVvariab_full_exec(symstructs, save_solutions, save_results)

%% Transient PhotoVoltage (TPV) - variable pulse intensity, various light intensities
%TPVvariab_full_exec(symstructs, false, false);

%% Transient PhotoVoltage (TPV) - variable pulse intensity, just dark & saving solutions
%TPVvariab_full_exec(ssol_i_eq, true, true);

%% Transient PhotoVoltage (TPV) - variable pulse intensity, just one sun
%TPVvariab_full_exec(ssol_i_light, false, false);

% test Impedance Spectroscopy with voltage step (ISstep)
% ISstep_struct = ISstep_full_exec(symstructs, deltaV, frozen_ions, save_solutions, save_results)

%% Impedance Spectroscopy with voltage step (ISstep) - various light intensities
%ISstep_full_exec(symstructs, 1e-3, 1, false, false, false);

%% Impedance Spectroscopy with voltage step (ISstep) - frozen ions & saving solutions
%ISstep_full_exec(ssol_i_light, 5e-4, 1, true, true, true);

%% Impedance Spectroscopy with voltage step (ISstep) - boundary conditions
% deltaV has to be very small for the solver not to go crazy
%ISstep_full_exec(ssol_i_light_BC2, 0.000015, 2, false, false, false);

%% Impedance Spectroscopy with voltage step (ISstep) - various deltaV
%ISstep_full_exec(ssol_i_light, [-1e-2, -1e-3, 1e-3, 1e-2], 1, false, false, false);

% test Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel)
% ISwave_results =  ISwave_full_exec_nonparallel(structs, startFreq, endFreq, Freq_points, deltaV, sequential, frozen_ions, do_graphics, save_solutions)

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC and various light intensities
ISwave_full_exec_nonparallel(symstructs, 1e2, 1e2, 1, 0.002, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, 1S and various frequencies
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 6, 0.002, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, dark and various frequencies
ISwave_full_exec_nonparallel(ssol_i_eq_SR, 1e8, 1e-2, 6, 0.002, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - SC and various light intensities
ISwave_full_exec_nonparallel(asymstructs, 1e2, 1e2, 1, 0.002, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, 1S and big oscillations
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e2, 1e2, 1, 0.2, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC, 1S and do not reach stability and do sequential simulations
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 100, 1, 3, 0.002, true, false, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - frozen ions & do graphics & saving solutions
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, false, true, true, true);

% test Impedance Spectroscopy with oscillating voltage (ISwave)
% ISwave_results = ISwave_full_exec(structs, startFreq, endFreq, Freq_points, deltaV, frozen_ions, do_graphics)

%% Impedance Spectroscopy with oscillating voltage (ISwave) - OC and various light intensities
ISwave_full_exec(symstructs, 1e2, 1e2, 1, 0.002, false, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - OC, dark and various frequencies
ISwave_full_exec(ssol_i_eq_SR, 1e8, 1e-2, 6, 0.002, false, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - frozen ions & do graphics
ISwave_full_exec(ssol_i_1S_SR, 1e8, 1e-2, 2, 0.002, true, true);

%% test Stark Effect - ElectroAbsorption (EA)
%true;