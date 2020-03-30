function soleq = equilibrate(varargin)
% Uses initial conditions defined in DF and runs to equilibrium
%
%% Input arguments
% VARARGIN{1,1} = PAR
% VARARGIN{1,2} = ELECTRONIC_ONLY
% ELECTRONIC_ONLY: 
% 0 = runs full equilibrate protocol
% 1 = skips ion equilibration
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
if length(varargin) == 1
    par = varargin{1,1};
    electronic_only = 0;
elseif length(varargin) == 2
    par = varargin{1,1};
    electronic_only = varargin{1,2};
else
    par = pc;
    electronic_only = 0;
end

tic;    % Start stopwatch
%% Initial arguments
% Setting sol.u = 0 enables a parameters structure to be read into
% DF but indicates that the initial conditions should be the
% analytical solutions
sol.u = 0;

% Store the original parameter set
par_origin = par;

% Start with zero SRH recombination
par.SRHset = 0;
% Radiative rec could initially be set to zero in addition if required
par.radset = 1;
% Start with no ionic carriers
par.N_ionic_species = 0;

%% General initial parameters
par.tmesh_type = 2;
par.tpoints = 10;

par.JV = 0;
par.Vapp = 0;
par.int1 = 0;
par.int2 = 0;
par.g1_fun_type = 'constant';
par.g2_fun_type = 'constant';
par.OC = 0;
par.tmesh_type = 2;
par.tmax = 1e-9;
par.t0 = par.tmax/1e4;
par.Rs = 0;
par.BC = 3;     % To enable boundary fluxes to be switched off

%% Switch off mobilities
par.mobset = 0;
par.mobseti = 0;

%% Initial solution with zero mobility
disp('Initial solution, zero mobility')
sol = df(sol, par);
disp('Complete')

% Switch on mobilities
par.BC = par_origin.BC;
par.mobset = 1;
par.radset = 1;

par.tmax = 1e-9;
par.t0 = par.tmax/1e6;

%% Solution with mobility switched on
disp('Solution with mobility switched on')
sol = df(sol, par);

par.tmax = 1e-3;
par.t0 = par.tmax/1e6;

sol = df(sol, par);

all_stable = verifyStabilization(sol.u, sol.t, 0.7);

% loop to check ions have reached stable config- if not accelerate ions by
% order of mag
j = 1;

while any(all_stable) == 0
    disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10^j)]);
    
    par.tmax = par.tmax*10;
    par.t0 = par.tmax/1e6;
    
    sol = df(sol, par);
    
    all_stable = verifyStabilization(sol.u, sol.t, 0.7);
end

soleq_nosrh = sol;

disp('Switching on interfacial recombination')
par.SRHset = 1;

par.tmax = 1e-6;
par.t0 = par.tmax/1e3;

soleq.el = df(soleq_nosrh, par);
disp('Complete')

if electronic_only == 0
    %% Equilibrium solutions with ion mobility switched on
    par.N_ionic_species = par_origin.N_ionic_species;
    
    % CReate temporary solution for appending initial conditions to
    sol = soleq_nosrh;
    
    % Write ionic initial conditions
    switch par.N_ionic_species
        case 1
            sol.u(:,:,4) = repmat(par.dev.Nani, length(sol.t), 1);
        case 2
            sol.u(:,:,4) = repmat(par.dev.Nani, length(sol.t), 1);
            sol.u(:,:,5) = repmat(par.dev.Nani, length(sol.t), 1);
    end
    
    % Start without SRH or series resistance
    par.SRHset = 0;
    par.Rs = 0;
%     
    disp('Closed circuit equilibrium with ions')
    
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
    
    par.mobset = 1;
    par.mobseti = 1;           % Ions are accelerated to reach equilibrium
    par.K_anion = rat_anion;
    par.K_cation = rat_cation;
    par.tmax = 1e-12;
    par.t0 = par.tmax/1e3;
      
    sol = df(sol, par);
    
    % Longer second solution
    par.calcJ = 0;
    par.tmax = 1e-3;
    par.t0 = par.tmax/1e3;
    
    sol = df(sol, par);
    
    all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    
    % loop to check ions have reached stable config- if not accelerate ions by
    % order of mag
    while any(all_stable) == 0
        disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10^j)]);
        
        par.tmax = par.tmax*10;
        par.t0 = par.tmax/1e6;
        
        sol = df(sol, par);
        
        all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    end
    
    % write solution and reset ion mobility
    soleq_i_nosrh = sol;
    soleq_i_nosrh.par.mobseti = 1;
    soleq_i_nosrh.par.K_anion = 1;
    soleq_i_nosrh.par.K_cation = 1;
    
    disp('Ion equilibrium solution complete')
    
    %% Ion equilibrium with surface recombination
    disp('Switching on SRH recombination')
    par.SRHset = 1;
    
    par.calcJ = 0;
    par.tmax = 1e-6;
    par.t0 = par.tmax/1e3;
    par.mobseti = 1;
    par.K_anion = 1;
    par.K_cation = 1;
    
    soleq.ion = df(soleq_i_nosrh, par);
    disp('Complete')
end

disp('EQUILIBRATION COMPLETE')
toc

end