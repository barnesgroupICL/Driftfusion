%% Intrinsic
par.standard = pc('input_files/spiro_mapi_tio2.csv');
soleq.standard = equilibrate(par.standard);

%% standard eppr = 5
par.standard_epp5 = pc('input_files/spiro_mapi_tio2_epp5.csv');
soleq.standard_epp5 = equilibrate(par.standard_epp5);

%% standard eppr = 50
par.standard_epp50 = pc('input_files/spiro_mapi_tio2_epp50.csv');
soleq.standard_epp50 = equilibrate(par.standard_epp50);

%% Open circuit
OC.standard_epp25 = lightonRs(soleq.standard.el, 1, -1e-2, 0, 1e6, 200);
OC.standard_epp5 = lightonRs(soleq.standard_epp5.el, 1, -1e-2, 0, 1e6, 200);
OC.standard_epp50 = lightonRs(soleq.standard_epp50.el, 1, -1e-2, 0, 1e6, 200);
OC.standard_i = lightonRs(soleq.standard.ion, 1, -1e-2, 1, 1e6, 200);

%% JVs
JV.standard_epp5 = doJV(soleq.standard_epp5.el, 1e-3, 100, 1, 0, 0, 1.2, 2);
JV.standard_epp25 = doJV(soleq.standard.el, 1e-3, 100, 1, 0, 0, 1.2, 2);
JV.standard_epp50 = doJV(soleq.standard_epp5.el, 1e-3, 100, 1, 0, 0, 1.2, 2);
JV.standard_i = doJV(soleq.standard.ion, 1e-6, 100, 1, 1, 0, 1.2, 2);