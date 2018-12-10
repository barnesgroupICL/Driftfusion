function soleq = equilibrate(varargin)
% Uses initial conditions defined in PINDRIFT and runs to equilibrium
if length(varargin) == 1
    par = varargin{1,1};
else
    par = pc;
end

tic;    % Start stopwatch
%% Initial arguments
% Setting sol.sol = 0 enables a parameters structure to be read into
% DF but indicates that the initial conditions should be the
% analytical solutions
sol.sol = 0;    

p_original = par;

% Store initial mobility of intrinsic layer- note all mobilities will be
% set to this value during the equilibration procedure.
sn_l = par.sn_l;
sn_r = par.sn_r;
sp_l = par.sp_l;
sp_r = par.sp_r;

BC = par.BC;

%% Start with low interfacial recombination coefficients
par.SRHset = 0;

% Raditative recombination could also be set to low values initially if required. 
% par.krad = 1e-20;
% par.kradetl = 1e-20;
% par.kradhtl = 1e-20;

%% General initial parameters
par.tmesh_type = 2;
par.tpoints = 40;

par.JV = 0;
par.Vapp = 0;
par.Int = 0;
par.pulseon = 0; 
par.OC = 0;
par.BC = BC;
par.tmesh_type = 2;
par.tmax = 1e-9;
par.t0 = par.tmax/1e4;

%% Switch off mobilities
par.mobset = 0;
par.mobseti= 0;

% Switch off extraction and recombination
par.sn_l = 0;
par.sn_r = 0;
par.sp_l = 0;
par.sp_r = 0;

%% Initial solution with zero mobility
disp('Initial solution, zero mobility')
sol = df(sol, par);
disp('Complete')

% Switch on mobilities
par.mobset = 1;

par.sn_l = p_original.sn_l;
par.sn_r = p_original.sn_r;
par.sp_l = p_original.sp_l;
par.sp_r = p_original.sp_r;

par.figson = 1;
par.tmax = 1e-9;
par.t0 = par.tmax/1e3;

%% Soluition with mobility switched on
disp('Solution with mobility switched on')
sol = df(sol, par);

% switch to desired mobilities
par.mue = p_original.mue;          % electron mobility
par.muh = p_original.muh;      % hole mobility

sol = df(sol, par);

par.tmax = 1e-3;
par.t0 = par.tmax/1e6;

sol = df(sol, par);

all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

% loop to check ions have reached stable config- if not accelerate ions by
% order of mag
j = 1;

while any(all_stable) == 0
    disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10^j)]);
    
    par.tmax = par.tmax*10;
    par.t0 = par.tmax/1e6;

    sol = df(sol, par);
    
    all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

end

soleq.eq = sol;

disp('Complete')

disp('Switching on interfacial recombination')
par.taun_inter = p_original.taun_inter;
par.taup_inter = p_original.taup_inter;
par.SRHset = 1;
par.dev = pc.builddev(par);

par.calcJ = 0;
par.tmax = 1e-6;
par.t0 = par.tmax/1e3;

soleq.eq_sr = df(soleq.eq, par);
disp('Complete')

%% Equilibrium solutions with ion mobility switched on

% Start with low recombination coefficients
par.SRHset = 0;

disp('Closed circuit equilibrium with ions')

par.mobseti = 1e4;           % Ions are accelerated to reach equilibrium
par.OC = 0;
par.tmax = 1e-9;
par.t0 = par.tmax/1e3; 

sol = df(soleq.eq, par);

% Much longer second step to ensure that ions have migrated
par.calcJ = 0;
par.tmax = 1e-2;
par.t0 = par.tmax/1e3;

sol = df(sol, par);

all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

% loop to check ions have reached stable config- if not accelerate ions by
% order of mag

while any(all_stable) == 0
    disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10^j)]);
    
    par.tmax = par.tmax*10;
    par.t0 = par.tmax/1e6;

    sol = df(sol, par);
    
    all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

end

% write solution and reset ion mobility
soleq.i = sol;
soleq.i.par.mobseti = 1;

disp('Ion equilibrium solution complete')

%% Ion equilibrium with surface recombination
disp('Switching on surface recombination')
par.taun_inter = p_original.taun_inter;
par.taup_inter = p_original.taup_inter;
par.SRHset = 1;
par.dev = pc.builddev(par);

par.calcJ = 0;
par.tmax = 1e-6;
par.t0 = par.tmax/1e3;

soleq.i_sr = df(soleq.i, par);
disp('Complete')

disp('EQUILIBRATION COMPLETE')
toc

end
