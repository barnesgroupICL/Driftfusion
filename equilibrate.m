function soleq = equilibrate(varargin)
% Uses analytical initial conditions and runs to equilibrium
% Note that tmax is consistently adjusted to appropriate values for to
% ensure there are numerous mesh points where large gradients in the time
% dimension are present
if length(varargin) == 1
    p = varargin{1,1};
else
    p = pc;
end

tic;    % Start stopwatch

%% Initial arguments
% Setting sol.sol = 0 enables a parameters structure to be read into
% pindrift but indicates that the initial conditions should be the
% analytical solutions
sol.sol = 0;    

p_original = p;

% Store initial mobility of intrinsic layer- note all mobilities will be
% set to this value during the equilibration procedure.

sn_l = p.sn_l;
sn_r = p.sn_r;
sp_l = p.sp_l;
sp_r = p.sp_r;

BC = p.BC;

%% Start with low interfacial recombination coefficients
p.taun_inter = [1e6, 1e6];       % [s] SRH time constant for electrons
p.taup_inter = [1e6, 1e6];      % [s] SRH time constant for holes

% Raditative recombination could also be set to low values initially if required. 
% p.krad = 1e-20;
% p.kradetl = 1e-20;
% p.kradhtl = 1e-20;

%% General initial parameters
p.tmesh_type = 2;
p.tpoints = 40;

p.JV = 0;
p.Vapp = 0;
p.Int = 0;
p.pulseon = 0; 
p.OC = 0;
p.BC = BC;
p.tmesh_type = 2;
p.tmax = 1e-9;
p.t0 = p.tmax/1e4;

%% Switch off mobilities
p.mobset = 0;
p.mobseti= 0;

% Switch off extraction and recombination
p.sn_l = 0;
p.sn_r = 0;
p.sp_l = 0;
p.sp_r = 0;

%% Initial solution with zero mobility
disp('Initial solution, zero mobility')
sol = pindrift(sol, p);
disp('Complete')

% Switch on mobilities
p.mobset = 1;

p.sn_l = p_original.sn_l;
p.sn_r = p_original.sn_r;
p.sp_l = p_original.sp_l;
p.sp_r = p_original.sp_r;

p.figson = 1;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;

%% Soluition with mobility switched on
disp('Solution with mobility switched on')
sol = pindrift(sol, p);

% switch to desired mobilities
p.mue = p_original.mue;          % electron mobility
p.muh = p_original.muh;      % hole mobility

sol = pindrift(sol, p);

p.tmax = 1e-3;
p.t0 = p.tmax/1e6;

sol = pindrift(sol, p);

all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

% loop to check ions have reached stable config- if not accelerate ions by
% order of mag
j = 1;

while any(all_stable) == 0
    disp(['increasing equilibration time, tmax = ', num2str(p.tmax*10^j)]);
    
    p.tmax = p.tmax*10;
    p.t0 = p.tmax/1e6;

    sol = pindrift(sol, p);
    
    all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

end

soleq.eq = sol;

disp('Complete')

disp('Switching on interfacial recombination')
p.taun_inter = p_original.taun_inter;
p.taup_inter = p_original.taup_inter;

p.calcJ = 0;
p.tmax = 1e-6;
p.t0 = p.tmax/1e3;

soleq.eq_sr = pindrift(soleq.eq, p);
disp('Complete')

%% Start with low recombination coefficients
p.taun_inter = [1e6, 1e6];       % [s] SRH time constant for electrons
p.taup_inter = [1e6, 1e6];      % [s] SRH time constant for holes

%% Equilibrium solutions with ion mobility switched on
%% Closed circuit conditions
disp('Closed circuit equilibrium with ions')

p.mobseti = 1;
p.OC = 0;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p.muion = [0,1e-6,0];           % Ions are accelerated to reach equilibrium
p.dev = pc.builddev(p);

sol = pindrift(soleq.eq, p);

% Much longer second step to ensure that ions have migrated
p.calcJ = 0;
p.tmax = 1e-2;
p.t0 = p.tmax/1e3;

sol = pindrift(sol, p);

all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

% loop to check ions have reached stable config- if not accelerate ions by
% order of mag
j = 1;

while any(all_stable) == 0
    disp(['Accelerating ions for equilibrium solution, mui = ', num2str(p.mui*10^j)]);
    
    % use mobseti as a multiplier- bit dangerous? But fine for now
    p.mobseti = p.mobseti*10;
    
    sol = pindrift(sol, p);
    
    all_stable = verifyStabilization(sol.sol, sol.t, 0.7);

end
sol.p.mobseti = 1;
% write solution and reset ion mobility
soleq.i = sol;

disp('Ion equilibrium solution complete')


%% Ion equilibrium with surface recombination
disp('Switching on surface recombination')
p.taun_inter = p_original.taun_inter;
p.taup_inter = p_original.taup_inter;

p.calcJ = 0;
p.tmax = 1e-6;
p.t0 = p.tmax/1e3;

soleq.i_sr = pindrift(soleq.i, p);
disp('Complete')

disp('EQUILIBRATION COMPLETE')
toc

end
