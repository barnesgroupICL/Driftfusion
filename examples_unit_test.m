% test pindrift and photophysics examples
% USAGE: runtests('unit_test')

% prepare solutions and indirectly test equilibrate and genIntStructs
% functions
[sol_eq, sol_i_eq, sol_i_eq_SR, ssol_i_eq, ssol_i_eq_SR, sol_i_1S_SR, ssol_i_1S_SR] = equilibrate_minimal();
symstructs = genIntStructs(ssol_i_eq_SR, 100, 1e-7, 4, true);
asymstructs = genIntStructs(sol_i_eq_SR, 100, 1e-7, 4, true);

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
% ISwave_struct =  ISwave_full_exec_nonparallel(structs, startFreq, endFreq, Freq_points, deltaV, sequential, frozen_ions, do_graphics, save_solutions)

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC and various light intensities
ISwave_full_exec_nonparallel(symstructs, 1e2, 1e2, 1, 0.002, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC and various frequencies and just one input solution at 1 sun
ISwave_full_exec_nonparallel(ssol_i_1S_SR, 1e8, 1e-2, 6, 0.002, 1, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - OC and various frequencies and just one input solution at 1 sun
ISwave_full_exec_nonparallel(ssol_i_eq_SR, 1e8, 1e-2, 6, 0.002, 1, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage and without parallelization (ISwave_nonparallel) - SC and various light intensities
ISwave_full_exec_nonparallel(asymstructs, 1e2, 1e2, 1, 0.002, 1, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - big oscillations
ISwave_full_exec(ssol_i_1S_SR, 500, 500, 1, 1.0, 1, false, false, false, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - boundary conditions
ISwave_full_exec(ssol_i_light_BC2, 1e9, 1e-1, 4, 0.002, 2, false, false, false, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - reach stability
ISwave_full_exec(ssol_i_light, 500, 500, 1, 0.002, 2, true, false, false, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - frozen ions & saving solutions
ISwave_full_exec(ssol_i_light, 1e9, 1e-1, 4, 0.002, 1, false, true, false, false, false, true, true);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - calculate ionic current
ISwave_full_exec(ssol_i_light, 500, 500, 1, 0.002, 1, false, false, true, false, false, false, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - parallelization
ISwave_full_exec(ssol_i_light, 1e9, 1e-1, 4, 0.002, 1, false, false, false, true, false, false, false);

%% Impedance Spectroscopy with oscillating voltage (ISwave) - do graphics switch
ISwave_full_exec(ssol_i_light, 500, 500, 1, 0.002, 1, false, false, false, false, true, false, false);

%% test Stark Effect - ElectroAbsorption (EA)

true;