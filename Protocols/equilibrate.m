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
% Switch off volumetric surface recombination check
par.vsr_check = 0;

%% General initial parameters
% Set applied bias to zero
par.V_fun_type = 'constant';
par.V_fun_arg(1) = 0;

% Set light intensities to zero
par.int1 = 0;
par.int2 = 0;
par.g1_fun_type = 'constant';
par.g2_fun_type = 'constant';

% Time mesh
par.tmesh_type = 2;
par.tpoints = 10;

% Series resistance
par.Rs = 0;

%% Switch electronic mobilities to 1
par.mu_n(par.mu_n > 0) = 1;
par.mu_p(par.mu_p > 0) = 1;
par.sn(par.sn > 0) = 1e-10;
par.sp(par.sp > 0) = 1e-10;
par.sn_l(par.sn_l > 0) = 1e5;     % electron surface recombination velocity left boundary
par.sn_r(par.sn_r > 0) = 1e5;     % electron extraction velocity right boundary
par.sp_l(par.sp_l > 0) = 1e5;     % hole extraction left boundary
par.sp_r(par.sp_r > 0) = 1e5;     % hole surface recombination velocity right boundary
par = refresh_device(par);

%% Switch off mobilities
par.mobset = 0;
par.mobseti = 0;

%% Initial solution with zero mobility
disp('Initial solution, zero mobility')
sol = df(sol, par);
disp('Complete')

% Switch on mobilities
par.mobset = 1;
par.radset = 1;
par.SRHset = 1;

% Characteristic diffusion time
t_diff = (par.dcum0(end)^2)/(2*par.kB*par.T*min(min(par.mu_n), min(par.mu_p)));
t_diff_norm = (par.dcum0(end)^2)/(2*par.kB*par.T);

par.tmax = 1000*t_diff_norm;
par.t0 = par.tmax/1e6;

%% Solution with mobility switched on
disp('Solution with mobility switched on')
sol = df(sol, par);

all_stable = verifyStabilization(sol.u, sol.t, 0.7);

% loop to check electrons have reached stable config- if not accelerate ions by
% order of mag
j = 1;

while any(all_stable) == 0
    disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10^j)]);

    par.tmax = 10*par.tmax;
    par.t0 = par.tmax/1e6;

    sol = df(sol, par);

    all_stable = verifyStabilization(sol.u, sol.t, 0.7);
end

%% reset mobilities
par.mu_n = par_origin.mu_n;
par.mu_p = par_origin.mu_p;
par.sn = par_origin.sn;
par.sp = par_origin.sp;
par.sn_l = par_origin.sn_l;    % electron surface recombination velocity left boundary
par.sn_r = par_origin.sn_r;     % electron extraction velocity right boundary
par.sp_l = par_origin.sp_l;     % hole extraction left boundary
par.sp_r = par_origin.sp_r;     % hole surface recombination velocity right boundary
par = refresh_device(par);

%% Run final solution with mobilties set back to defined
par.tmax = 10*t_diff;
par.t0 = par.tmax/1e6;

sol = df(sol, par);

soleq.el = sol;
% Manually check final section of solution for VSR self-consitency
sol_ic = extract_IC(soleq.el, [soleq.el.t(end)*0.7, soleq.el.t(end)]);
compare_rec_flux(sol_ic, par.RelTol_vsr, par.AbsTol_vsr, 0);
% Switch VSR check on for future use
soleq.el.par.vsr_check = 1;

disp('Electronic carrier equilibration complete')

if electronic_only == 0 && par_origin.N_ionic_species > 0
    %% Equilibrium solutions with ion mobility switched on
    par.N_ionic_species = par_origin.N_ionic_species;

    % Create temporary solution for appending initial conditions to
    sol = soleq.el;

    % Start without SRH or series resistance
    %par.SRHset = 0;
    par.Rs = 0;

    disp('Closed circuit equilibrium with ions')

    % Take ratio of electron and ion mobilities in the active layer
    par.mu_n(par.mu_n > 0) = 1;
    par.mu_p(par.mu_p > 0) = 1;
    par.sn(par.sn > 0) = 1e-10;
    par.sp(par.sp > 0) = 1e-10;
    par.mu_c(par.mu_c > 0) = 1;
    par.mu_a(par.mu_a > 0) = 1;
    par.sn_l(par.sn_l > 0) = 1e5;     % electron surface recombination velocity left boundary
    par.sn_r(par.sn_r > 0) = 1e5;     % electron extraction velocity right boundary
    par.sp_l(par.sp_l > 0) = 1e5;     % hole extraction left boundary
    par.sp_r(par.sp_r > 0) = 1e5;     % hole surface recombination velocity right boundary
    
    par = refresh_device(par);
    
    par.mobset = 1;
    par.mobseti = 1;           % Ions are accelerated to reach equilibrium
    par.tmax = 100*t_diff_norm;
    par.t0 = par.tmax/1e6;

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
    
    %% Reset mobilities and boundary coefficients
    par.mu_n = par_origin.mu_n;
    par.mu_p = par_origin.mu_p;
    par.mu_c = par_origin.mu_c;
    par.mu_a = par_origin.mu_a;
    par.sn = par_origin.sn;
    par.sp = par_origin.sp;
    par.sn_l = par_origin.sn_l;    % electron surface recombination velocity left boundary
    par.sn_r = par_origin.sn_r;     % electron extraction velocity right boundary
    par.sp_l = par_origin.sp_l;     % hole extraction left boundary
    par.sp_r = par_origin.sp_r;     % hole surface recombination velocity right boundary
    par = refresh_device(par);
    
    %% Final solution with parameters set back to defined ones
    par.tmax = 10*t_diff;
    par.t0 = par.tmax/1e6;
    sol = df(sol, par);
    
    % write solution
    soleq.ion = sol;
    % Manually check solution for VSR self-consitency
    sol_ic = extract_IC(soleq.ion, [soleq.ion.t(end)*0.7, soleq.ion.t(end)]);
    compare_rec_flux(sol_ic, par.RelTol_vsr, par.AbsTol_vsr, 0);
    % Reset switches
    soleq.ion.par.vsr_check = 1;
    soleq.ion.par.mobseti = 1;
    soleq.ion.par.K_a = 1;
    soleq.ion.par.K_c = 1;

    disp('Ionic carrier equilibration complete')
end

disp('EQUILIBRATION COMPLETE')
toc

end
