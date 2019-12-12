initialise_df

%% Read in base parameters
par_da3l = pc('Input_files/spiro_mapi_tio2_3lda.csv');

%% Get equilibrium solutions
soleq_da3l = equilibrate(par_da3l);

%% Do a 20 mV jump forward and then record the decay
% Input arguments: jumptoV(sol_ini, Vjump, tdwell, mobseti, Int, stabilise, accelerate)
sol_jump_da3l = jumptoV(soleq_da3l.ion, 200e-3, 10, 1, 0, 1, 0);

%% Extract active layer properties from params
al = par_da3l.active_layer;
pcum0 = par_da3l.pcum0;
pmid = pcum0(al) + round((pcum0(al+1)-pcum0(al))/2);

%% Fit the decay of the cation current
Jepp10 = dfana.calcJ(sol_jump_da3l);

p1 = find(sol_jump_da3l.t < sol_jump_da3l.t(end)/10);
p1_fit = p1(end);
p2 = find(sol_jump_da3l.t < sol_jump_da3l.t(end)/2);
p2_fit = p2(end);

fit_DF = fit(sol_jump_da3l.t(p1_fit:p2_fit)',log(Jepp10.c(p1_fit:p2_fit, pmid)),'poly1');
figure;plot(fit_DF,sol_jump_da3l.t,log(Jepp10.c(:,pmid)));

%% Read out the time constants DF = Driftfusion
tau_DF = -1/fit_DF.p1;

%% Plot Electrostatic profile for Driftfusion for comparison with DA
dfplot.Vx(sol_jump_da3l, [0,0.1, 0.2, 0.4, 1.6]);
hold on

%% Get time constants from Depletion Approximation using 20 mV step forward in Voltage DA = Depletion Approximation
% depletion_approx_modelX_Vjump(par, tmax, Vstart, Vjump, tpoints, tarr, figson)
[rhomat, Fmat, Phimat, Q, tout] = depletion_approx_modelX_Vjump(par_da3l, 10, 0, 200e-3, 1000, [0,0.1, 0.2, 0.4, 1.6], 1);

%% Fit the curves
p1 = find(tout < tout(end)/10);
p1_fit = p1(end);
p2 = find(tout < tout(end)/4);
p2_fit = p2(end);

deltaQ = abs(Q)- min(abs(Q));
fit_DA = fit(tout(p1_fit:p2_fit), log(deltaQ(p1_fit:p2_fit)'), 'poly1');
tau_DA = -1/(fit_DA.p1);
figure;plot(fit_DA, tout, log(deltaQ));
xlabel('Time [s]')
ylabel('log(Q)')

