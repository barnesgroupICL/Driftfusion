%% Read in params
par_1l = pc('Input_files/1_layer_intrinsic.csv');
par_pn_epp3 = pc('Input_files/2_layer_pn_epp3.csv');
par_pn_epp12 = pc('Input_files/2_layer_pn_epp12.csv');
par_pn_epp48 = pc('Input_files/2_layer_pn_epp48.csv');

%% Get equilibrium
soleq_1l = equilibrate(par_1l);
soleq_pn_epp3 = equilibrate(par_pn_epp3);
soleq_pn_epp12 = equilibrate(par_pn_epp12);
soleq_pn_epp48 = equilibrate(par_pn_epp48);

%% Illuminate OC
OC_1l = lightonRs(soleq_1l.el, 1, 1, 0, 1e8, 200);
OC_pn_epp3 = lightonRs(soleq_pn_epp3.el, 1, 1, 0, 1e8, 200);
OC_pn_epp12 = lightonRs(soleq_pn_epp12.el, 1, 1, 0, 1e8, 200);
OC_pn_epp48 = lightonRs(soleq_pn_epp48.el, 1, 1, 0, 1e8, 200);
OC_pn_epp12_ion = lightonRs(soleq_pn_epp12.ion, 1, 10000, 1, 1e8, 200);

%% Calculate recombination and kpfo
G_1l = dfana.calcG(OC_1l);
G_pn_epp3 = dfana.calcG(OC_pn_epp3);
G_pn_epp12 = dfana.calcG(OC_pn_epp12);
G_pn_epp48 = dfana.calcG(OC_pn_epp48);
G_pn_epp12_ion = dfana.calcG(OC_pn_epp12_ion);

R_1l = dfana.calcR(OC_1l);
R_pn_epp3 = dfana.calcR(OC_pn_epp3);
R_pn_epp12 = dfana.calcR(OC_pn_epp12);
R_pn_epp48 = dfana.calcR(OC_pn_epp48);
R_pn_epp12_ion = dfana.calcR(OC_pn_epp12_ion);

[kpfo_n_OC_1l, n_integrated_OC_1l] = dfana.calckpfo_n(OC_1l);
[kpfo_n_OC_pn_epp3, n_integrated_OC_pn_epp3] = dfana.calckpfo_n(OC_pn_epp3);
[kpfo_n_OC_pn_epp12, n_integrated_OC_pn_epp12] = dfana.calckpfo_n(OC_pn_epp12);
[kpfo_n_OC_pn_epp48, n_integrated_OC_pn_epp48] = dfana.calckpfo_n(OC_pn_epp48);
[kpfo_n_OC_pn_epp12_ion, n_integrated_OC_pn_epp12_ion] = dfana.calckpfo_n(OC_pn_epp12_ion);

Voc_1l = abs(dfana.calcVQFL(OC_1l));
Voc_1l = Voc_1l(end);
Voc_pn_epp3 = abs(dfana.calcVQFL(OC_pn_epp3));
Voc_pn_epp3 = Voc_pn_epp3(end);
Voc_pn_epp12 = abs(dfana.calcVQFL(OC_pn_epp12));
Voc_pn_epp12 = Voc_pn_epp12(end);
Voc_pn_epp48 = abs(dfana.calcVQFL(OC_pn_epp48));
Voc_pn_epp48 = Voc_pn_epp48(end);
Voc_pn_epp12_ion = abs(dfana.calcVQFL(OC_pn_epp12_ion));
Voc_pn_epp12_ion = Voc_pn_epp12_ion(end);

%% Plots equilibrium
dfplot.ELnpx(soleq_1l.el);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.ELnpx(soleq_pn_epp3.el);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.ELnpx(soleq_pn_epp12.el);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.ELnpx(soleq_pn_epp48.el);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
legend('','intrinsic, n', 'intrinsic, p','pn, epp=3, n',...
    'pn, epp=3 p','pn epp=12, n', 'pn epp=12, p', 'pn epp=48, n', 'pn epp=48, p')
hold off

%% Energy level and carrier OC
dfplot.ELnpx(OC_1l);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.ELnpx(OC_pn_epp3);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.ELnpx(OC_pn_epp12);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.ELnpx(OC_pn_epp48);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
legend('','intrinsic, n', 'intrinsic, p','pn, epp=3, n',...
    'pn, epp=3 p','pn, epp=12, n', 'pn, epp=12, p', 'pn, epp=48, n', 'pn, epp=48, p')

%% EL,np, OC ions vs no ions
dfplot.ELx(OC_pn_epp12);
subplot(3, 1, 1)
hold on
subplot(3, 1, 2)
hold on
subplot(3, 1, 3)
hold on
dfplot.ELx(OC_pn_epp12_ion);
subplot(3, 1, 1)
hold off
subplot(3, 1, 2)
hold off
legend('','','pn, epp=12, NSD=0 cm-3, n', 'pn, epp=12, NSD=0 cm-3, p',...
    'pn, epp=12, NSD=1e19 cm-3, n', 'pn, epp=12, NSD=1e19 cm-3, p')
subplot(3, 1, 3)
hold off

%% rhox Fx Vx
dfplot.rhoxFxVx(OC_1l);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.rhoxFxVx(OC_pn_epp3);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.rhoxFxVx(OC_pn_epp12);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
dfplot.rhoxFxVx(OC_pn_epp48);
subplot(2, 1, 1)
hold on
subplot(2, 1, 2)
hold on
legend('','intrinsic', 'pn, epp=3', 'pn, epp=12', 'pn, epp=48', 'pn, epp=12, NSD=1e19 cm-3')

%% rhox
dfplot.rhox(OC_1l);
hold on
dfplot.rhox(OC_pn_epp3);
hold on
dfplot.rhox(OC_pn_epp12);
hold on
dfplot.rhox(OC_pn_epp48);
hold on
dfplot.rhox(OC_pn_epp12_ion);
hold on
legend('','intrinsic', 'pn, epp=3', 'pn, epp=12', 'pn, epp=48', 'pn, epp=12, NSD=1e19 cm-3')

%% Vx
dfplot.Vx(OC_1l);
hold on
dfplot.Vx(OC_pn_epp3);
hold on
dfplot.Vx(OC_pn_epp12);
hold on
dfplot.Vx(OC_pn_epp48);
hold on
dfplot.Vx(OC_pn_epp12_ion);
hold on
legend('','intrinsic', 'pn, epp=3', 'pn, epp=12', 'pn, epp=48', 'pn, epp=12, NSD=1e19 cm-3')

