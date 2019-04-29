function JV = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
tic
disp('Current voltage scan')

% A procedure for running JV scans using DF
%% Input arguments
% sol_ini   	= an initial solution - must be at the same Vapp as the
% starting voltage for the scan- currently set to Vapp = 0 V
% JVscan_rate   = self-explanatory
% JVscan_pnts   = no. of points in time mesh
% mobseti       = determines whether ion mobility is on or off
% Vstart        = scan start voltage
% Vend          = scan end voltage
% option        1 = dark only, 2 = light only, 3 = dark & light

%[JV.dk.f, JV.dk.r, JV.ill.f, JV.ill.r] = deal(0);

% Read parameters structure into structure par
par = sol_ini.par;

%%Initial settings
par.Ana = 1;
par.figson = 1;
par.Int = 0;
par.pulseon = 0;
par.OC = 0;
par.mobseti = mobseti;

%% JV settings
par.JV = 1;
par.Vstart = Vstart;
par.Vend = Vend;
par.calcJ = 0;
par.JVscan_pnts = JVscan_pnts;
par.tmax = abs(par.Vend- par.Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
par.t0 = 0;
par.tmesh_type = 1;
par.tpoints = par.JVscan_pnts;

if option ==1 || option ==3
        
        %% Dark forward scan
        
        disp('Dark forward scan...')
        JV.dk.f = df(sol_ini, par);
        JV.dk.f.par.JV = 0;
        disp('Complete.')
        
        %% Dark reverse scan
        
        disp('Dark reverse scan...')
        par.Vstart = Vend;
        par.Vend = Vstart;
        par.JV = 1;
        
        JV.dk.r = df(JV.dk.f, par);
        JV.dk.r.par.JV = 0;
        disp('Complete.')
        
end

if option ==2 || option ==3
        
        %% 1 Sun quasi equilibrium solution
        
        disp('1 Sun quasi-equilibrium solution')
        par.JV = 0;
        par.mobseti = 0;          % Switch ion mobility off for illumination step
        par.Int = Intensity;
        
        % Log time mesh
        par.tmesh_type = 2;
        par.tmax = 1e-3;
        par.t0 = par.tmax*1e-6;
        
        sol_i_1S = df(sol_ini, par);
        disp('Complete.')
        
        %% Light forward
        
        disp('Light forward scan...')
        par.mobseti = mobseti;
        
        %% JV settings
        par.JV = 1;
        par.Vstart = Vstart;
        par.Vend = Vend;
        par.calcJ = 0;
        par.tmax = abs(par.Vend- par.Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
        par.t0 = 0;
        par.tmesh_type = 1;
        par.tpoints = par.JVscan_pnts;
        
        JV.ill.f = df(sol_i_1S, par);
        JV.ill.f.par.JV = 0;
        disp('Complete.')
        
        %% Light reverse
      
        disp('Light reverse scan...')
        par.Vstart = Vend;
        par.Vend = Vstart;
        
        JV.ill.r = df(JV.ill.f, par);
        JV.ill.r.par.JV = 0;
        disp('Complete.')
        
        figure(11)
        hold off
        disp('JV scan complete.')
        
        plotJV(JV, option)
        JV.stats = JVstats(JV);
end
toc
end