function JVsol = doJV(sol_ini, JVscan_rate, JVscan_pnts, Intensity, mobseti, Vstart, Vend, option)
% A procedure for running JV scans using DF
% Input arguments
% sol_ini   	= an initial solution - must be at the same Vapp as the
% starting voltage for the scan- currently set to Vapp = 0 V
% JVscan_rate   = self-explanatory
% JVscan_pnts   = no. of points in time mesh
% mobseti       = determines whether ion mobility is on or off
% Vstart        = scan start voltage
% Vend          = scan end voltage
% option        1 = dark only, 2 = light only, 3 = dark & light
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
tic
disp('Current voltage scan')

% Read parameters structure into structure par
par = sol_ini.par;

%%Initial settings
par.int1 = 0;
par.mobseti = mobseti;

%% JV settings
par.tmesh_type = 1;
par.t0 = 0;
par.tmax = abs(Vend - Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
par.tpoints = JVscan_pnts;

par.V_fun_type = 'sweep';
par.V_fun_arg(1) = Vstart;
par.V_fun_arg(2) = Vend;
par.V_fun_arg(3) = par.tmax;

if option ==1 || option ==3

    %% Dark forward scan

    disp('Dark forward scan...')
    JVsol.dk.f = df(sol_ini, par);
    disp('Complete.')

    %% Dark reverse scan
    disp('Dark reverse scan...')

    % Sweep settings
    par.V_fun_arg(1) = getVend(JVsol.dk.f);
    par.V_fun_arg(2) = Vstart;
    par.V_fun_arg(3) = par.tmax;

    JVsol.dk.r = df(JVsol.dk.f, par);
    disp('Complete.')

end

if option ==2 || option ==3

    if Intensity == 0
        error('Intensity cannot be zero- use option 1 for dark scan instead')
    end

    par = sol_ini.par;      % reset parameters
    %% 1 Sun quasi equilibrium solution
    disp('Illuminated quasi-equilibrium solution')
    par.tmesh_type = 1;
    par.tmax = 1e-3;
    par.t0 = 0;
    par.tpoints = 10;

    par.mobseti = 0;          % Switch ion mobility off for illumination step
    % Set up light sweep
    par.g1_fun_type = 'sweep';
    % COEFF = [Amplitude_initial, Amplitude_final, tmax]
    par.g1_fun_arg = [0, Intensity, par.tmax];

    par.V_fun_type = 'constant';
    par.V_fun_arg(1) = Vstart;

    sol_i_1S = df(sol_ini, par);
    disp('Complete.')

    %% Light forward
    par.g1_fun_type = 'constant';
    % COEFF = [Amplitude_initial, Amplitude_final, tmax]
    par.g1_fun_arg = Intensity;
    par.int1 = Intensity;

    disp('Illuminated forward scan...')
    par.mobseti = mobseti;

    %% JV settings
    par.tmesh_type = 1;
    par.t0 = 0;
    par.tmax = abs(Vend- Vstart)/JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
    par.tpoints = JVscan_pnts;

    % Sweep settings
    par.V_fun_type = 'sweep';
    par.V_fun_arg(1) = Vstart;
    par.V_fun_arg(2) = Vend;
    par.V_fun_arg(3) = par.tmax;

    JVsol.ill.f = df(sol_i_1S, par);
    disp('Complete.')

    %% Light reverse
    disp('Illuminated reverse scan...')
    par.V_fun_arg(1) = getVend(JVsol.ill.f);
    par.V_fun_arg(2) = Vstart;

    JVsol.ill.r = df(JVsol.ill.f, par);
    disp('Complete.')

    disp('JV scan complete.')

    JVsol.stats = dfana.JVstats(JVsol);
end
toc
end
