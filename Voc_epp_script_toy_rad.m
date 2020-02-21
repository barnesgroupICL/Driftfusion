%% Intrinsic
par.symmetric = pc('input_files/3_layer_test_symmetric.csv');
soleq.symmetric = equilibrate(par.symmetric);

%% symmetric eppr = 5
par.symmetric_epp5 = pc('input_files/3_layer_test_symmetric_epp5.csv');
soleq.symmetric_epp5 = equilibrate(par.symmetric_epp5);

%% symmetric eppr = 50
par.symmetric_epp50 = pc('input_files/3_layer_test_symmetric_epp50.csv');
soleq.symmetric_epp50 = equilibrate(par.symmetric_epp50);

%% Open circuit
OC.symmetric_epp25 = lightonRs(soleq.symmetric.el, 1, -1e-2, 0, 1e6, 200);
OC.symmetric_epp5 = lightonRs(soleq.symmetric_epp5.el, 1, -1e-2, 0, 1e6, 200);
OC.symmetric_epp50 = lightonRs(soleq.symmetric_epp50.el, 1, -1e-2, 0, 1e6, 200);
OC.symmetric_i = lightonRs(soleq.symmetric.ion, 1, -1e-2, 1, 1e6, 200);

%% JVs
JV.symmetric_epp5 = doJV(soleq.symmetric_epp5.el, 1e-3, 121, 1, 0, 0, 1.2, 2);
JV.symmetric_epp25 = doJV(soleq.symmetric.el, 1e-3, 121, 1, 0, 0, 1.2, 2);
JV.symmetric_epp50 = doJV(soleq.symmetric_epp5.el, 1e-3, 121, 1, 0, 0, 1.2, 2);
JV.symmetric_i = doJV(soleq.symmetric.ion, 1e-6, 121, 1, 1, 0, 1.2, 2);