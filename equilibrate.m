function soleq = equilibrate(varargin)
% Uses analytical initial conditions and runs to equilibrium
% Note that tmax is consistently adjusted to appropriate values for to
% ensure there are numerous mesh points where large gradients in the time
% dimension are present
if length(varargin) == 1
    p = varargin{1,1};
else
    disp('Building parameters object...')
    p = pc;
    disp('Complete')
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

%% Start with SRH switched off
p.SRHset = 0;

% Raditative recombination could also be set to low values initially if required. 
% p.krad = 1e-20;
% p.kradetl = 1e-20;
% p.kradhtl = 1e-20;

%% General initial parameters
p.tmesh_type = 2;
p.tpoints = 40;

p.Ana = 1;
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
p.mobseti = 0;

% Switch off extraction and recombination
p.sn_l = 0;
p.sn_r = 0;
p.sp_l = 0;
p.sp_r = 0;

%% Initial solution with zero mobility
disp('Initial solution, zero mobility')
sol = pindrift(sol, p);
disp('Complete')

% Switch on mobilities - easy numbers first
p.mobset = 1;
p.mobseti = 0;

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

p.Ana = 1;
p.calcJ = 0;
p.tmax = 1e-2;
p.t0 = p.tmax/1e10;

soleq.eq = pindrift(sol, p);
disp('Complete')

%% Equilibrium solutions with ion mobility switched on
%% Closed circuit conditions
disp('Closed circuit equilibrium with ions')
p.mobseti = 1;

p.OC = 0;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p.muion = [0, 1e-6, 0];           % Ions are accelerated to reach equilibrium
p.dev = pc.builddev(p);

sol = pindrift(soleq.eq, p);

% Much longer second step to ensure that ions have migrated
p.calcJ = 0;
p.tmax = 1e-2;
p.t0 = p.tmax/1e3;

soleq.i = pindrift(sol, p);
% Switch to default ion mobility
soleq.i.p.muion = p_original.muion;
soleq.i.p.dev = pc.builddev(p_original);

disp('Complete')

%% Ion equilibrium with surface recombination
disp('Switching on surface recombination')
p.SRHset = 1;

p.calcJ = 0;
p.tmax = 1e-6;
p.t0 = p.tmax/1e3;

soleq.i_sr = pindrift(soleq.i, p);
disp('Complete')
% Switch to default ion mobility
soleq.i_sr.p.muion = p_original.muion;
soleq.i_sr.p.dev = soleq.i.p.dev;

% % Switch off SR
% p.taun = [1e6, 1e6, 1e6];
% p.taup = [1e6, 1e6, 1e6];
% 
% %% Symmetricise closed circuit condition
% disp('Symmetricise equilibriumion solution')
% symsol = symmetricize(soleq.i);
% disp('Complete')
% 
% p.OC = 1;
% p.tmax = 1e-9;
% p.t0 = p.tmax/1e3;
% 
% %% Switch off mobilities
% p.mue = [0, 0, 0];          % electron mobility
% p.muh = [0, 0, 0];      % hole mobility
% 
% %% OC condition with ions at equilbirium
% disp('Open circuit equilibrium with ions')
% ssol = pindrift(symsol, p);
% 
% p.tmax = 1e-9;
% p.t0 = p.tmax/1e3;
% 
% % Switch on mobilities
% p.mue = p_original.mue;          % electron mobility
% p.muh = p_original.muh;      % hole mobility
% p.mui = 0;
% 
% ssol = pindrift(ssol, p);
% 
% % Switch on ion mobility to ensure equilibrium has been reached
% p.tmax = 1e-9;
% p.t0 = p.tmax/1e3;
% p.mui = 1e-6;
% 
% ssol = pindrift(ssol, p);
% 
% p.tmax = 1e-2;
% p.t0 = p.tmax/1e3;
% 
% ssol_i_eq = pindrift(ssol, p);
% 
% disp('Complete')
% 
% %% Ions, OC Surface recombination
% p.taun = p_original.taun;
% p.taup = p_original.taup;
% 
% p.tmax = 1e-3;
% p.t0 = p.tmax/1e6;
% 
% ssol_i_eq_SR = pindrift(ssol_i_eq , p);
% 
% %% 1 Sun quasi equilibrium
% disp('1 Sun quasi equilibrium')
% tmax = 1e-3;
% p.t0 = p.tmax/1e6;
% p.Int = 1;
% 
% ssol_i_1S_SR = pindrift(ssol_i_eq_SR, p);
% 
% disp('Complete')

%}
%}
disp('EQUILIBRATION COMPLETE')
toc

end
