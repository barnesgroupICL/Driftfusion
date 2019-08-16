%% A script to perform a transients of the transient photovoltage
% measurement
% close all
% clear all

% Read-in params
par_tio2 = pc('input_files/spiro_mapi_tio2.csv');

% Find equilibrium solution
soleq_tio2 = equilibrate(par_tio2);

% Perform open circiut voltage transient
% lighton_Rs(sol_ini, int1, stab_time, mobseti, Rs, pnts)
sol_OC = lighton_Rs(soleq_tio2.ion, 1, 10, 1, 1e6, 400);

% Create array of Ntr linearly spaced times
Ntr = 6;    % Number of voltage transients
tarr = 0:sol_OC.par.tmax/Ntr:sol_OC.par.tmax;

% Create temporary solution that solutions at specific times will be
% written into
sol_temp = sol_OC;

%% Loop to get transients
figure(101)
hold on

for i = 1:Ntr
    % Find the approximate time point
    pp = find(sol_OC.t <= tarr(i));
    pp = pp(end);
    
    % Write solution into final position of temporary solution (this is
    % what will be read in
    sol_temp.u = sol_OC.u(pp,:,:);
    % Reduce the time step for greater accuracy
    sol_temp.par.MaxStepFactor = 0.1;
    
    % Freeze ions during pulse, 1 us pulse, with 49 us decay
    %sol_pulse] = dolightpulse(sol_ini, pulse_int, tmax, tpoints, duty, mobseti, log_timemesh)
    sol_VTROTTR(i) = do_light_pulse(sol_temp, 1, 50e-6, 400, 2, 0, 1);
    
    % Extract Voc vs t and subtract baseline
    Voc(i,:) = dfana.Voct(sol_VTROTTR(i));
    deltaV(i,:) = Voc(i,:) -  Voc(i,1);
    
    % Plot the outputs
    plot(sol_VTROTTR(i).t, deltaV(i,:))
           
end
xlabel('Time [s]')
ylabel('DeltaV [V]')
xlim([0, 1e-5])
hold off

% Plot the slow Voc transient
dfplot.Voct(sol_OC);

% Plot generation profile as a function of x_ihalf and t
x_ihalf = getvarihalf(sol_VTROTTR(1).x);

%% Plot the generation rate as a function of x and t
dfplot.gxt(sol_VTROTTR(1))