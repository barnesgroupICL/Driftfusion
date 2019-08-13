function [sol_Voc, Voc] = findVoc(sol_ini, Int, mobseti, x0, x1, tol)
% FINDVOC finds the open cicuit voltage of a device using a
% Newton-Raphson iteration

%% Instructions
% Before running the code, run a JV scan to get the approximate Voc value.
% Then use values either side of the Voc to obtain good intial guesses for
% x0 and x1.

%% Input arguments
% SOL_INI- initial state of the device
% MOBSETI - ion mobility switch
% X0 - lower initial guess
% X1 - upper initial guess
% TOL - current tolerence in Acm-2 - the value that the current should drop to below
% before the procedure stops running

par = sol_ini.par;

%% Initital settings
par.OC = 0;
par.calcJ = 0;
par.pulseon = 0;
JVscan_rate = 1;        %Vs-1

%% Save Surface rec coefficients
taun = par.taun_inter;       % [s] SRH time constant for electrons
taup = par.taup_inter;

%% 1 Sun quasi equilibrium solution
disp('1 Sun quasi-equilibrium solution')
par.mobseti = 0;      % switch off ion mobility for illumination step
par.JV = 0;
par.int1 = 1;

% Log time mesh
par.tmesh_type = 2;
par.tmax = 1e-3;
par.t0 = par.tmax*1e-6;

sol = df(sol_ini, par);
disp('Complete.')

%% Run JV to new potential
par.mobseti = mobseti;
par.Ana = 1;
par.JV = 1;
par.Vstart = par.Vapp;
par.Vend = x0;
par.calcJ = 0;
par.tmax = abs(par.Vend- par.Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
par.t0 = 0;
par.tmesh_type = 1;
par.tpoints = par.JVscan_pnts;
par.pulseon = 0;
par.int1 = Int;
      
sol = df(sol, par);         

%% Stabilised solution
disp('Initial stabilised solution')
par.JV = 0;
par.Vapp = x0;
par.tmesh_type = 2;
par.tmax = 1e-2;
par.t0 = par.tmax/1e4;
%par = mobsetfun(1, 0, par); 

sol = df(sol, par);

disp('Complete')

fx0 = sol.Jtotr(end, end);

xrun = x0;
yrun = fx0;

fx1 = 1;%                 % Set any initial vale > 0.01

figure(111)

while abs(fx1) > tol
    
    i = 1;
    
    disp(['Newton-Raphson iteration', num2str(i)])
 
    %% Scan to new potential x1
    disp('Scan to new potential x1')
    
    par.JV = 1;
    par.Vstart = x0;
    par.Vend = x1;
    par.tmax = abs(par.Vend- par.Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
    par.t0 = 0;
    par.tmesh_type = 1;
    par.tpoints = par.JVscan_pnts;

    sol = df(sol, par);         
    
    disp('Complete')
    
    %% Stabilised solution
    
    disp('Stabilised solution')
    
    par.JV = 0;
    par.Vapp = x1;
    par.tmesh_type = 2;
    par.tmax = 1e-2;
    par.t0 = par.tmax/1e4;
    %par = mobsetfun(1, 0, par); 

    sol = df(sol, par);
    
    disp('Complete')
    fx1 = sol.Jtotr(end, end);
% 
%     par.tmesh_type = 2;
%     par.tmax = 1e-9;
%     par.t0 = par.tmax/1e3;
%     par.Vapp = x1;
%     par.pulseon = 0;
%     par.int1 = 0;
% 
%     sol = df(sol_eq, par);         
%     
%     par.BC = 1;
%     par.tmax = 1e-6;
%     par.t0 = par.tmax/1e4;
%     par = mobsetfun(1, 0, par); 
% 
%     sol = df(sol, par); 
%     
%     % Switch on surface rec
%     par.taun_etl = 1e6;       % [s] SRH time constant for electrons
%     par.taup_etl = par.taun_etl;      % [s] SRH time constant for holes
%     par.taun_htl = par.taun_etl; 
%     par.taup_htl = par.taun_etl; 
%     
%     sol = df(sol, par);
%     
%     par.int1 = Int;
%     par.tmax = 1e-3;
%     par.t0 = par.tmax/1e4;
% 
%     sol = df(sol, par);
% 
%     par = mobsetfun(1, mobi, par); 
% 
%     par.tmax = 1e0;
%     par.t0 = par.tmax/1e4;
% 
%     sol = df(sol, par);
% 
%     par.tmax = 1e-6;
%     par.t0 = par.tmax/1e4;
%     par.calcJ = 2;
% 
%     sol_V1 = df(sol, par);
% 
%     fx1 = sol_V1.Jn(end);
    
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
    
end %endwhile

sol_Voc = sol;

Voc = x1;

end