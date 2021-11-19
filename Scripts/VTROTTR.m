% A script to perform a transients of the transient photovoltage
% measurement
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
%% Read-in params
par_tio2 = pc('input_files/spiro_mapi_tio2.csv');

%% Find equilibrium solution
soleq_tio2 = equilibrate(par_tio2);

%% Perform open circiut voltage transient
% lightonRs(sol_ini, int1, stab_time, mobseti, Rs, pnts)
sol_OC = lightonRs(soleq_tio2.ion, 1, 10, 1, 1e6, 400);

%% Create array of Ntr linearly spaced times
Ntr = 6;    % Number of voltage transients
tarr = 0:sol_OC.par.tmax/Ntr:sol_OC.par.tmax;

%% Create temporary solution that solutions at specific times will be
% written into
sol_temp = sol_OC;

%% Loop to get transients
figure(101)
hold on

for i = 1:Ntr
    %% Find the approximate time point
    time_point = find(sol_OC.t <= tarr(i));
    time_point = time_point(end);
    
    %% Write solution into final position of temporary solution (this is
    % what will be read in
    sol_temp.u = sol_OC.u(time_point,:,:);
    
    %% Reduce the time step for greater accuracy
    sol_temp.par.MaxStepFactor = 0.1;
    
    %% Freeze ions during pulse, 1 us pulse, with 49 us decay
    % sol_pulse = doLightPulse(sol_ini, pulse_int, tmax, tpoints, duty, mobseti, log_timemesh)
    sol_VTROTTR(i) = doLightPulse(sol_temp, 1, 50e-6, 400, 2, 0, 1);
    
    % Extract Voc vs t and subtract baseline
    Voc(i,:) = dfana.calcVQFL(sol_VTROTTR(i));
    deltaV(i,:) = Voc(i,:) -  Voc(i,1);
    
    % Plot the outputs
    plot(sol_VTROTTR(i).t, deltaV(i,:))
           
end
xlabel('Time [s]')
ylabel('DeltaV [V]')
xlim([0, 1e-5])
hold off

%% Plot the slow Voc transient
dfplot.Voct(sol_OC);

%% Plot the generation rate as a function of x and t
dfplot.gxt(sol_VTROTTR(1))