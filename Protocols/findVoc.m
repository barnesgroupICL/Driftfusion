function [sol_Voc, Voc] = findVoc(sol_ini, Int, mobseti, x0, x1, tol)
% FINDVOC finds the open cicuit voltage of a device using a
% Newton-Raphson iteration

%% Instructions
% Before running the code, it is abvisable to run a JV scan to get the approximate Voc value.
% Then use values either side of the Voc to obtain good intial guesses for
% x0 and x1.

%% Input arguments
% SOL_INI- initial state of the device
% MOBSETI - ion mobility switch
% X0 - lower initial guess
% X1 - upper initial guess
% TOL - current tolerence in Acm-2 - the value that the current should drop to below
% before the procedure stops running
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
par = sol_ini.par;

%% Initital settings
if mobseti == 1
    % Take ratio of electron and ion mobilities in the active layer
    rat_anion = par.mue(par.active_layer)/par.muani(par.active_layer);
    rat_cation = par.mue(par.active_layer)/par.mucat(par.active_layer);
    
    % If the ratio is infinity (ion mobility set to zero) then set the ratio to
    % zero instead
    if isnan(rat_anion) || isinf(rat_anion)
        rat_anion = 0;
    end
    
    if isnan(rat_cation) || isinf(rat_cation)
        rat_cation = 0;
    end
    
    par.mobseti = 1;           % Ions are accelerated to reach equilibrium
    par.K_anion = rat_anion;
    par.K_cation = rat_cation;
end

%% 1 Sun quasi equilibrium solution
disp('1 Sun quasi-equilibrium solution')
par.mobseti = 0;      % switch off ion mobility for illumination step
par.int1 = 1;

% Log time mesh
par.tmesh_type = 2;
par.tmax = 1e-3;
par.t0 = par.tmax*1e-6;

sol = df(sol_ini, par);
disp('Complete.')

%% Run JV to new potential
par.mobseti = mobseti;
par.V_fun_type = 'sweep';
par.V_fun_arg(1) = par.Vapp;
par.V_fun_arg(2) = x0;
par.V_fun_arg(3) = par.tmax;
par.t0 = 0;
par.tmesh_type = 1;
par.tpoints = 100;
par.int1 = Int;

disp('Voltage sweep to initial guess')
sol = df(sol, par);
disp('Complete.')

%% Stabilised solution
par.V_fun_type = 'constant';
par.V_fun_arg(1) = x0;
par.tmesh_type = 2;
par.tmax = 1e-3;
par.t0 = par.tmax/1e4;

disp('Initial stabilised solution')
sol = df(sol, par);
disp('Complete')

J = dfana.calcJ(sol);
fx0 = J.tot(end, end);

xrun = x0;
yrun = fx0;

fx1 = 1;%                 % Set any initial vale > 0.01

figure(111)

while abs(fx1) > tol
    
    i = 1;
    
    disp(['Newton-Raphson iteration', num2str(i)])
    
    %% Scan to new potential x1
    par.V_fun_type = 'sweep';
    par.V_fun_arg(1) = x0;
    par.V_fun_arg(2) = x1;
    par.V_fun_arg(3) = par.tmax;
    par.tmesh_type = 1;
    par.tmax = 1e-3;
    par.t0 = 0;
    
    disp('Scan to new potential x1')
    sol = df(sol, par);
    disp('Complete')
    
    %% Stabilised solution
    disp('Stabilised solution')
    
    par.V_fun_type = 'constant';
    par.V_fun_arg(1) = x1;
    par.tmesh_type = 2;
    par.tmax = 1e-3;
    par.t0 = par.tmax/1e4;
    
    sol = df(sol, par);
    
    disp('Complete')
    
    J = dfana.calcJ(sol);
    fx1 = J.tot(end, end);
    
    % approximation to the gradient
    fgrad = (fx1 - fx0)/(x1-x0)
    
    xplot = 0.6:0.01:1.6;
    yplot = fgrad*xplot + fx0- fgrad*x0;
    
    xrun = [xrun, x1];
    yrun = [yrun, fx1];
    
    plot(xrun, yrun, 'o', xplot, yplot);
    xlim([0.6, 1.2]);
    %ylim([-2e-6, 2e-6])
    xlabel('x');
    ylabel('fx');
    grid on;
    
    drawnow update
    % New initial guesses
    x0 = x1
    
    if fgrad ~= 0
        
        x1 = x1 - (fx1/fgrad);
        
    else
        
        x1 = x1 + 0.1;
        
    end
    
    fx0 = fx1
    %sol_V0 = sol_V1;
    
    i = i+1;
    
end

sol_Voc = sol;
sol_Voc.par.K_anion = 1;
sol_Voc.par.K_cation = 1;

Voc = x1;

end