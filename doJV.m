function JV = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
disp('Current voltage scan')

% A procedure for running JV scans using pindrift
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

% Read parameters structure into structure P
p = sol_ini.p;

%%Initial settings
p.Ana = 1;
p.figson = 1;
p.Int = 0;
p.pulseon = 0;
p.mobseti = mobseti;

%% JV settings
p.JV = 1;
p.Vstart = Vstart;
p.Vend = Vend;
p.calcJ = 0;
p.JVscan_pnts = JVscan_pnts;
p.tmax = abs(p.Vend- p.Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
p.t0 = 0;
p.tmesh_type = 1;
p.tpoints = p.JVscan_pnts;

if option ==1 || option ==3
        
        %% Dark forward scan
        
        disp('Dark forward scan...')
        JV.dk.f = pindrift(sol_ini, p);
        
        disp('Complete.')
        
        %% Dark reverse scan
        
        disp('Dark reverse scan...')
        p.Vstart = Vend;
        p.Vend = Vstart;
        p.JV = 1;
        
        JV.dk.r = pindrift(JV.dk.f, p);
        
        disp('Complete.')
        
end

if option ==2 || option ==3
        
        %% 1 Sun quasi equilibrium solution
        
        disp('1 Sun quasi-equilibrium solution')
        p.JV = 0;
        p.mobseti = 0;          % Switch ion mobility off for illumination step
        p.Int = Intensity;
        
        % Log time mesh
        p.tmesh_type = 2;
        p.tmax = 1e-3;
        p.t0 = p.tmax*1e-6;
        
        sol_i_1S = pindrift(sol_ini, p);
        disp('Complete.')
        
        %% Light forward
        
        disp('Light forward scan...')
        p.mobseti = mobseti;
        
        %% JV settings
        p.JV = 1;
        p.Vstart = Vstart;
        p.Vend = Vend;
        p.calcJ = 0;
        p.tmax = abs(p.Vend- p.Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
        p.t0 = 0;
        p.tmesh_type = 1;
        p.tpoints = p.JVscan_pnts;
        
        JV.ill.f = pindrift(sol_i_1S, p);
        
        disp('Complete.')
        
        %% Light reverse
      
        disp('Light reverse scan...')
        p.Vstart = Vend;
        p.Vend = Vstart;
        
        JV.ill.r = pindrift(JV.ill.f, p);
        disp('Complete.')
        
        figure(11)
        hold off
        disp('JV scan complete.')
        
        plotJV(JV, option)
        
end

JV.stats = JVstats(JV);

end

