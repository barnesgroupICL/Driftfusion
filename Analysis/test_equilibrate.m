function soleq = test_equilibrate(varargin)
% Uses initial conditions defined in DF and runs to equilibrium
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
% Radiative rec could initially be set to zero in addition if required
par.radset = 1;

%% General initial parameters
par.tmesh_type = 2;
par.tpoints = 40;

par.JV = 0;
par.Vapp = 0;
par.int1 = 0;
par.int2 = 0;
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

%% Initial solution with zero mobility
disp('Initial solution, zero mobility')
soleq = df(sol, par);
disp('Complete')
toc
end
