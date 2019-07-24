function soleq = equilibrate(varargin)
% Uses initial conditions defined in PINDRIFT and runs to equilibrium
if length(varargin) == 1
    par = varargin{1,1};
else
    par = pc;
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

% Raditative recombination could also be set to low values initially if required. 
% par.krad = 1e-20;
% par.kradetl = 1e-20;
% par.kradhtl = 1e-20;

%% General initial parameters
par.tmesh_type = 2;
par.tpoints = 10;

par.JV = 0;
par.Vapp = 0;
par.Int = 0;
par.pulseon = 0; 
par.OC = 0;
par.tmesh_type = 2;
par.tmax = 1e-9;
par.t0 = par.tmax/1e4;
par.Rs = 0;
par.Ana = 0;

%% Switch off mobilities
par.mobset = 0;
par.mobseti = 0;

% Switch off extraction and recombination
par.sn_l = 0;
par.sn_r = 0;
par.sp_l = 0;
par.sp_r = 0;

if par.DP3layerIC
 
    disp('Initial solution, zero mobility')
    sol = df(sol, par);
    disp('Complete')
    
    % Take ratio of electron and ion mobilities in the active layer
rat_anion = par.mue(par.active_layer)/par.muion(par.active_layer);
rat_cation = par.mue(par.active_layer)/par.mucat(par.active_layer);

% If the ratio is infinity (ion mobility set to zero) then set the ratio to
% zero instead
if isnan(rat_anion) || isinf(rat_anion)
    rat_anion = 0;
end

if isnan(rat_cation) || isinf(rat_cation)
    rat_cation = 0;
end

% Switch on mobilities and extraction
par.sn_l = par_origin.sn_l;
par.sn_r = par_origin.sn_r;
par.sp_l = par_origin.sp_l;
par.sp_r = par_origin.sp_r;
par.mobset = 1;
par.mobseti = 1;           % Ions are accelerated to reach equilibrium
par.K_anion = rat_anion;
par.K_cation = rat_cation;
par.tmax = 1e-3;
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


if par_origin.Rs ~= 0
    disp('Switching on series resistance')
    
    par.Rs = par_origin.Rs;
    par.tmax = 1e-6;
    par.t0 = 1e-12;
    
    soleq_nosrh = df(sol, par);
    disp('Complete')
end

sol = sol;

disp('Switching on interfacial recombination')
par.SRHset = 1;

par.tmax = 1e-6;
par.t0 = par.tmax/1e3;

soleq.ion = df(sol, par);
disp('Complete')
    
else


%% Initial solution with zero mobility
disp('Initial solution, zero mobility')
sol = df(sol, par);
disp('Complete')

% Switch on mobilities
par.mobset = 1;

par.sn_l = par_origin.sn_l;
par.sn_r = par_origin.sn_r;
par.sp_l = par_origin.sp_l;
par.sp_r = par_origin.sp_r;

par.tmax = 1e-9;
par.t0 = par.tmax/1e3;

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

if par_origin.Rs ~= 0
    disp('Switching on series resistance')
    
    par.Rs = par_origin.Rs;
    par.tmax = 1e-6;
    par.t0 = 1e-12;
    
    soleq_nosrh = df(sol, par);
    disp('Complete')
end

soleq_nosrh = sol;

disp('Switching on interfacial recombination')
par.SRHset = 1;

par.tmax = 1e-6;
par.t0 = par.tmax/1e3;

soleq.no_ion = df(soleq_nosrh, par);
disp('Complete')

%% Equilibrium solutions with ion mobility switched on

% Start without SRH or series resistance
par.SRHset = 0;
par.Rs = 0;

disp('Closed circuit equilibrium with ions')

% Take ratio of electron and ion mobilities in the active layer
rat_anion = par.mue(par.active_layer)/par.muion(par.active_layer);
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
par.tmax = 1e-12;
par.t0 = par.tmax/1e3; 

sol = df(soleq_nosrh, par);

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

disp('Switching on series resistance')

par.Rs = par_origin.Rs;  
par.tmax = 1e-6;
par.t0 = 1e-12;

sol = df(sol, par);

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

dfplot.ELx(soleq.ion);
dfplot.Jx(soleq.ion)

disp('EQUILIBRATION COMPLETE')
toc

end