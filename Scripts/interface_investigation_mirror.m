%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%% Start code
% Initialise the system
initialise_df

%% Create a par_mirrorameters object for Spiro/MAPI/TiO2 by including a filepath to the 
% appropriate .csv as the arugment to the par_mirrorameters class PC
par_mirror_mirror = pc('Input_files/3_layer_test_mirror.csv');

%% Good transport, no rec
par_mirror_ideal = par_mirror_mirror;
par_mirror_ideal.sn(2) = 1e-20;
par_mirror_ideal.sp(2) = 1e-20;
par_mirror_ideal.sn(4) = 1e-20;
par_mirror_ideal.sp(4) = 1e-20;
par_mirror_ideal.mue(2) = 1e3;
par_mirror_ideal.muh(2) = 1e3;
par_mirror_ideal.mue(4) = 1e3;
par_mirror_ideal.muh(4) = 1e3;
par_mirror_ideal = refresh_device(par_mirror_ideal);

sol_mirroreq_ideal = equilibrate(par_mirror_ideal);

% %sol_mirror_CV = doCV(sol_mirror_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
% sol_mirror_CV_ideal = doCV(sol_mirroreq_ideal.el, 1, 0, 1, 0, 100e-3, 1, 281);
%% With ions
sol_mirror_CV_ideal_ion = doCV(sol_mirroreq_ideal.ion, 1, 0, -1, 0, 100e-3, 1, 281);
%% With ions
sol_mirror_CV_ideal_ion_dk = doCV(sol_mirroreq_ideal.ion, 0, 0, -1, 0, 100e-3, 1, 281);

%% Poor transport, no rec
par_mirror_trans = par_mirror_mirror;
par_mirror_trans.sn(2) = 1e-20;
par_mirror_trans.sp(2) = 1e-20;
par_mirror_trans.sn(4) = 1e-20;
par_mirror_trans.sp(4) = 1e-20;
par_mirror_trans.mue(2) = 1e-5;
par_mirror_trans.muh(2) = 1e-4;
par_mirror_trans.mue(4) = 1e-3;
par_mirror_trans.muh(4) = 1e-6;
par_mirror_trans = refresh_device(par_mirror_trans);

sol_mirroreq_trans = equilibrate(par_mirror_trans);

% %sol_mirror_CV = doCV(sol_mirror_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
% sol_mirror_CV_trans = doCV(sol_mirroreq_trans.el, 1, 0, 1, 0, 100e-3, 1, 281);
% %%
% sol_mirror_CV_trans_dk = doCV(sol_mirroreq_trans.el, 0, 0, 1, 0, 100e-3, 1, 281);
%% With ions
sol_mirror_CV_trans_ion = doCV(sol_mirroreq_trans.ion, 1, 0, -1, 0, 100e-3, 1, 281);
%%
sol_mirror_CV_trans_ion_dk = doCV(sol_mirroreq_trans.ion, 0, 0, -1, 0, 100e-3, 1, 281);

%% High mob, High rec
par_mirror_rec = par_mirror_mirror;
par_mirror_rec.r_constant = 1e23;  % Uniform rec rate within interfacial regions
par_mirror_rec.sn(2) = 1e-20;      
par_mirror_rec.sp(2) = 1e-20;
par_mirror_rec.sn(4) = 1e-20;
par_mirror_rec.sp(4) = 1e-20;
par_mirror_rec.mue(2) = 1e3;
par_mirror_rec.muh(2) = 1e3;
par_mirror_rec.mue(4) = 1e3;
par_mirror_rec.muh(4) = 1e3;
par_mirror_rec = refresh_device(par_mirror_rec);

sol_mirroreq_rec = equilibrate(par_mirror_rec);

%sol_mirror_CV = doCV(sol_mirror_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
% %% Illuminated JV - electronic carriers only
% sol_mirror_CV_rec = doCV(sol_mirroreq_rec.el, 1, 0, 1, 0, 100e-3, 1, 281);
% %% Dark JV  - electronic carriers only
% sol_mirror_CV_rec_dk = doCV(sol_mirroreq_rec.el, 0, 0, 1, 0, 100e-3, 1, 281);
%% Illuminated JV - ionic and electronic
sol_mirror_CV_rec_ion = doCV(sol_mirroreq_rec.ion, 1, 0, -1, 0, 100e-3, 1, 281);
%% Dark JV - ionic and electronic
sol_mirror_CV_rec_ion_dk = doCV(sol_mirroreq_rec.ion, 0, 0, -1, 0, 100e-3, 1, 281);

%% Low mob, High rec - also asymmetric mob for testing
par_mirror_lomob_hirec = par_mirror_mirror;
par_mirror_lomob_hirec.r_constant = 1e23;  % Uniform rec rate within interfacial regions
par_mirror_lomob_hirec.sn(2) = 1e-20;      
par_mirror_lomob_hirec.sp(2) = 1e-20;
par_mirror_lomob_hirec.sn(4) = 1e-20;
par_mirror_lomob_hirec.sp(4) = 1e-20;
par_mirror_lomob_hirec.mue(2) = 1e-2;
par_mirror_lomob_hirec.muh(2) = 1e-4;
par_mirror_lomob_hirec.mue(4) = 1e-2;
par_mirror_lomob_hirec.muh(4) = 1e-4; 
par_mirror_lomob_hirec = refresh_device(par_mirror_lomob_hirec);

sol_mirroreq_lomob_hirec = equilibrate(par_mirror_lomob_hirec);

%sol_mirror_CV = doCV(sol_mirror_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
%% Illuminated JV - electronic carriers only
% sol_mirror_CV_lomob_hirec = doCV(sol_mirroreq_lomob_hirec.el, 1, 0, 1, 0, 100e-3, 1, 281);
%% Dark JV - electronic carriers only
% sol_mirror_CV_lomob_hirec_dk = doCV(sol_mirroreq_lomob_hirec.el, 0, 0, 1, 0, 100e-3, 1, 281); 
%% Illuminated JV - ionic and electronic
sol_mirror_CV_lomob_hirec_ion = doCV(sol_mirroreq_lomob_hirec.ion, 1, 0, -1, 0, 100e-3, 1, 281);
%% Dark JV - ionic and electronic
sol_mirror_CV_lomob_hirec_ion_dk = doCV(sol_mirroreq_lomob_hirec.ion, 0, 0, -1, 0, 100e-3, 1, 281);

% 
% %% Compar_mirrorison with analytical sol_mirrorutions
% compare_carrier_interfaces_7(sol_mirror_CV_ideal, [0, 4, 8, 12]);
% %%
% compare_carrier_interfaces_7(sol_mirror_CV_trans, [0, 4, 8, 12]);
% %%
% compare_carrier_interfaces_7(sol_mirror_CV_trans_dk, [0, 4, 8, 12]);
% %%
% compare_carrier_interfaces_7(sol_mirror_CV_rec, [0, 4, 8, 12]);
% %%
% compare_carrier_interfaces_7(sol_mirror_CV_rec_dk, [0, 4, 8, 12]);
% %%
% compare_carrier_interfaces_7(sol_mirror_CV_lomob_hirec, [0, 4, 8, 12]);

%% Compar_mirrorison with analytical sol_mirrorutions
%% DARK
compare_carrier_interfaces_X(sol_mirror_CV_ideal_ion_dk, [0, 4, 8, 12]);
%%
compare_carrier_interfaces_X(sol_mirror_CV_trans_ion_dk, [0, 4, 8, 12]);
%%
compare_carrier_interfaces_X(sol_mirror_CV_rec_ion_dk, [0, 4, 8, 12]);
%%
compare_carrier_interfaces_X(sol_mirror_CV_lomob_hirec_ion, [0, 4, 8, 12]);

%% ILLUMINATED
%%
compare_carrier_interfaces_X(sol_mirror_CV_ideal_ion, [0, 4, 8, 12]);
%%
compare_carrier_interfaces_X(sol_mirror_CV_trans_ion, [0, 4, 8, 12]);
%%
compare_carrier_interfaces_X(sol_mirror_CV_rec_ion, [0, 4, 8, 12]);
%%
compare_carrier_interfaces_X(sol_mirror_CV_lomob_hirec_ion_dk, [0, 4, 8, 12]);
%% Plots
% dfplot.JtotVapp(sol_mirror_CV_rec,0)
% hold off
% ylim([-30e-3,10e-3])
% dfplot.JtotVapp(sol_mirror_CV_ideal,0)
% hold on
% dfplot.JtotVapp(sol_mirror_CV_trans,0)
% hold on
% dfplot.JtotVapp(sol_mirror_CV_rec,0)
% hold off
% ylim([-30e-3,10e-3])

% dfplot.npx(sol_mirror_CV_ideal,5)
% hold on
% dfplot.npx(sol_mirror_CV_trans,5)
% hold on
% dfplot.npx(sol_mirror_CV_rec,5)
% hold off

