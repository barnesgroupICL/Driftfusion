function [JV_dk_f, JV_dk_r, JV_1S_f, JV_1S_r] = doJV(sol_ini, JVscan_rate, mui, Vstart, Vend)

% A procedure for running JV scans using pindrift
%% Input arguments
% sol_ini   	= an initial solution - must be at the same Vapp as the
% starting voltage for the scan- currently set to Vapp = 0 V
% JVscan_rate   = self-explanatory
% mui           = ion mobility
% Vstart        = Start voltage
% Vend          = End voltage

% Read parameters structure into structure P
p = sol_ini.p;

%%Initial settings
p.Ana = 1;
p.figson = 1;
p.Int = 0;
p.pulseon = 0;
p.mui = mui;

%% JV settings
p.JV = 1;
p.Vstart = Vstart;
p.Vend = Vend;
p.calcJ = 0;
p.tmax = abs(p.Vend- p.Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
p.t0 = 0;
p.tmesh_type = 1;
p.tpoints = p.JVscan_pnts;

%% Dark forward scan

disp('Dark forward scan...')
JV_dk_f = pindrift(sol_ini, p);

figure(11)
xlim([0, 1.3]);
ylim([-20, 20]);
hold on

disp('Complete.')

%% Dark reverse scan

disp('Dark reverse scan...')
p.Vstart = Vend;
p.Vend = Vstart;
p.JV = 1;

JV_dk_r = pindrift(JV_dk_f, p);

disp('Complete.')
%% 1 Sun quasi equilibrium solution

disp('1 Sun quasi-equilibrium solution')
p.JV = 0;
p.mui = 0;          % Switch ion mobility off for illumination step
p.Int = 1;

% Log time mesh
p.tmesh_type = 2;
p.tmax = 1e-3;
p.t0 = p.tmax*1e-6;

sol_i_1S = pindrift(sol_ini, p);
disp('Complete.')

%% Light forward

disp('Light forward scan...')
p.mui = mui;

%% JV settings
p.JV = 1;
p.Vstart = Vstart;
p.Vend = Vend;
p.calcJ = 0;
p.tmax = abs(p.Vend- p.Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
p.t0 = 0;
p.tmesh_type = 1;
p.tpoints = p.JVscan_pnts;

JV_1S_f = pindrift(sol_i_1S, p);

disp('Complete.')
%% Light reverse

disp('Light forward scan...')
p.Vstart = Vend;
p.Vend = Vstart;

JV_1S_r = pindrift(JV_1S_f, p);
disp('Complete.')

figure(11)
hold off
disp('JV scan complete.')
end

