function [sol_Voc, Voc] = findVoc(sol_ini, mobseti, x0, x1, tol)
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

p = sol_ini.p;

%% Initital settings
p.OC = 0;
p.calcJ = 0;
p.pulseon = 0;
Suns = 1;               % Intensity
JVscan_rate = 1;        %Vs-1

%% Save Surface rec coefficients
taun = p.taun;       % [s] SRH time constant for electrons
taup = p.taup;

%% 1 Sun quasi equilibrium solution
disp('1 Sun quasi-equilibrium solution')
p.mobseti = 0;      % switch off ion mobility for illumination step
p.JV = 0;
p.Int = 1;

% Log time mesh
p.tmesh_type = 2;
p.tmax = 1e-3;
p.t0 = p.tmax*1e-6;

sol = pindrift(sol_ini, p);
disp('Complete.')

%% Run JV to new potential
p.mobseti = mobseti;
p.Ana = 1;
p.JV = 1;
p.Vstart = p.Vapp;
p.Vend = x0;
p.calcJ = 0;
p.tmax = abs(p.Vend- p.Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
p.t0 = 0;
p.tmesh_type = 1;
p.tpoints = p.JVscan_pnts;
p.pulseon = 0;
p.Int = Suns;
      
sol = pindrift(sol, p);         

%% Stabilised solution
disp('Initial stabilised solution')
p.JV = 0;
p.Vapp = x0;
p.tmesh_type = 2;
p.tmax = 1e-2;
p.t0 = p.tmax/1e4;
%p = mobsetfun(1, 0, p); 

sol = pindrift(sol, p);

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
    
    p.JV = 1;
    p.Vstart = x0;
    p.Vend = x1;
    p.tmax = abs(p.Vend- p.Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
    p.t0 = 0;
    p.tmesh_type = 1;
    p.tpoints = p.JVscan_pnts;

    sol = pindrift(sol, p);         
    
    disp('Complete')
    
    %% Stabilised solution
    
    disp('Stabilised solution')
    
    p.JV = 0;
    p.Vapp = x1;
    p.tmesh_type = 2;
    p.tmax = 1e-2;
    p.t0 = p.tmax/1e4;
    %p = mobsetfun(1, 0, p); 

    sol = pindrift(sol, p);
    
    disp('Complete')
    fx1 = sol.Jtotr(end, end);
% 
%     p.tmesh_type = 2;
%     p.tmax = 1e-9;
%     p.t0 = p.tmax/1e3;
%     p.Vapp = x1;
%     p.pulseon = 0;
%     p.Int = 0;
% 
%     sol = pindrift(sol_eq, p);         
%     
%     p.BC = 1;
%     p.tmax = 1e-6;
%     p.t0 = p.tmax/1e4;
%     p = mobsetfun(1, 0, p); 
% 
%     sol = pindrift(sol, p); 
%     
%     % Switch on surface rec
%     p.taun_etl = 1e6;       % [s] SRH time constant for electrons
%     p.taup_etl = p.taun_etl;      % [s] SRH time constant for holes
%     p.taun_htl = p.taun_etl; 
%     p.taup_htl = p.taun_etl; 
%     
%     sol = pindrift(sol, p);
%     
%     p.Int = Suns;
%     p.tmax = 1e-3;
%     p.t0 = p.tmax/1e4;
% 
%     sol = pindrift(sol, p);
% 
%     p = mobsetfun(1, mobi, p); 
% 
%     p.tmax = 1e0;
%     p.t0 = p.tmax/1e4;
% 
%     sol = pindrift(sol, p);
% 
%     p.tmax = 1e-6;
%     p.t0 = p.tmax/1e4;
%     p.calcJ = 2;
% 
%     sol_V1 = pindrift(sol, p);
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