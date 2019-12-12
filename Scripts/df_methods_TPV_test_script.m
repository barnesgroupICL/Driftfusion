%% Get params
par_tpv = pc('Input_files/TPV_test.csv');

%% Get equilibrium
soleq_tpv = equilibrate(par_tpv);

% Decrease maximum timestep for better time accuracy (note slows solution
% considerably)
soleq_tpv.el.par.MaxStepFactor = 0.01;

%% Do the TPV
% Note in this instance due to the use of a single slab of material and
% perfectly blocking contacts, application of series resistance is not
% required
% doTPV(sol_ini, bias_int, stab_time, mobseti, Rs, pulse_int, tmax, tpoints, duty)
tmax = 20e-6;
duty = 5; 

[sol_TPV_0p1sun, sol_ill_0p1sun] = doTPV(soleq_tpv.el, 0.1, 10, 0, 0, 0.02, tmax, 1000, duty);
[sol_TPV_1sun, sol_ill_1sun] = doTPV(soleq_tpv.el, 1, 10, 0, 0, 0.2, tmax, 1000, duty);
[sol_TPV_10sun, sol_ill_10sun] = doTPV(soleq_tpv.el, 10, 10, 0, 0, 2, tmax, 1000, duty);

%% get Delta Voc
Voc_0p1sun = dfana.Voct(sol_TPV_0p1sun);
deltaVoc_0p1sun = Voc_0p1sun-Voc_0p1sun(1);

Voc_1sun = dfana.Voct(sol_TPV_1sun);
deltaVoc_1sun = Voc_1sun-Voc_1sun(1);

Voc_10sun = dfana.Voct(sol_TPV_10sun);
deltaVoc_10sun = Voc_10sun-Voc_10sun(1);

%% Shift time
t_TPV_0p1sun = sol_TPV_0p1sun.t - (tmax*(1e-2*duty));
t_TPV_1sun = sol_TPV_1sun.t - (tmax*(1e-2*duty));
t_TPV_10sun = sol_TPV_10sun.t - (tmax*(1e-2*duty));

%% Plot TPV
figure(600)
plot(t_TPV_0p1sun , deltaVoc_0p1sun, t_TPV_1sun , deltaVoc_1sun, t_TPV_10sun , deltaVoc_10sun)
xlabel('Time [s]')
ylabel('\Delta Voc')
hold on

%% Semilog plot TPV
figure(601)
semilogy(t_TPV_0p1sun , deltaVoc_0p1sun, t_TPV_1sun , deltaVoc_1sun, t_TPV_10sun , deltaVoc_10sun)
xlabel('Time [s]')
ylabel('\Delta Voc')
hold on

%% Obtain zero dimensional solution
TPV0Dstruct_1sun = TPV_0D(par_tpv, 0.1, 0.02, tmax, duty);
TPV0Dstruct_1sun = TPV_0D(par_tpv, 1, 0.2, tmax, duty);
TPV0Dstruct_1sun = TPV_0D(par_tpv, 10, 2, tmax, duty);

figure(600)
legend('DF, 0.1 sun', 'DF, 1 sun', 'DF, 10 sun', 'Analytical')
hold off
figure(601)
legend('DF, 0.1 sun', 'DF, 1 sun', 'DF, 10 sun', 'Analytical')
hold off
