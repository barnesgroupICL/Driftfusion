%% symmetric_srh eppr = 5
par.symmetric_srh_epp5 = pc('input_files/3_layer_test_symmetric_srh_epp5.csv');
soleq.symmetric_srh_epp5 = equilibrate(par.symmetric_srh_epp5);

%% symmetric, epp = 25
par.symmetric_srh_epp25  = pc('input_files/3_layer_test_symmetric_srh_epp25.csv');
soleq.symmetric_srh_epp25  = equilibrate(par.symmetric_srh_epp25 );

%% symmetric_srh eppr = 50
par.symmetric_srh_epp50 = pc('input_files/3_layer_test_symmetric_srh_epp50.csv');
soleq.symmetric_srh_epp50 = equilibrate(par.symmetric_srh_epp50);

%% Open circuit
OC.symmetric_srh_epp5 = lightonRs(soleq.symmetric_srh_epp5.el, 1, -1e-2, 0, 1e8, 200);
Voct.epp5 = dfana.QFLs(OC.symmetric_srh_epp5);
Voct.epp5 = abs(Voct.epp5(end));
[OC_Newton.symmetric_srh_epp5, Voc.symmetric_srh_epp5] =...
    findVoc(soleq.symmetric_srh_epp5.el, 1, 0, Voct.epp5*0.9, Voct.epp5*1.1, 1e-12);

OC.symmetric_srh_epp25 = lightonRs(soleq.symmetric_srh_epp25.el, 1, -1e-2, 0, 1e8, 200);
Voct.epp25 = dfana.QFLs(OC.symmetric_srh_epp25);
Voct.epp25 = abs(Voct.epp25(end));
[OC_Newton.symmetric_srh_epp25, Voc.symmetric_srh_epp25] =...
    findVoc(soleq.symmetric_srh_epp25.el, 1, 0, Voct.epp25*0.9, Voct.epp25*1.1, 1e-12);

OC.symmetric_srh_epp50 = lightonRs(soleq.symmetric_srh_epp50.el, 1, -1e-2, 0, 1e8, 200);
Voct.epp50 = dfana.QFLs(OC.symmetric_srh_epp50);
Voct.epp50 = abs(Voct.epp50(end));
[OC_Newton.symmetric_srh_epp50, Voc.symmetric_srh_epp50] =...
    findVoc(soleq.symmetric_srh_epp50.el, 1, 0, Voct.epp50*0.9, Voct.epp50*1.1, 1e-12);
%%
OC.symmetric_srh_epp25_i = lightonRs(soleq.symmetric_srh_epp25.ion, 1, -1e-2, 1, 1e8, 200);
Voct.epp25_i = dfana.QFLs(OC.symmetric_srh_epp25_i);
Voct.epp25_i = abs(Voct.epp25_i(end));
[OC_Newton.symmetric_srh_epp25_i, Voc.symmetric_srh_epp25_i] =...
    findVoc(soleq.symmetric_srh_epp25.ion, 1, 1, Voct.epp25_i*0.9, Voct.epp25_i*1.1, 1e-12);

%% JVs
JV.symmetric_srh_epp5 = doJV(soleq.symmetric_srh_epp5.el, 1e-3, 121, 1, 0, 0, 1.2, 2);
JV.symmetric_srh_epp25 = doJV(soleq.symmetric_srh_epp25.el, 1e-3, 121, 1, 0, 0, 1.2, 2);
JV.symmetric_srh_epp50 = doJV(soleq.symmetric_srh_epp50.el, 1e-3, 121, 1, 0, 0, 1.2, 2);
JV.symmetric_srh_i = doJV(soleq.symmetric_srh_epp25.ion, 1e-6, 121, 1, 1, 0, 1.2, 2);