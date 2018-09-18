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

%% Start with low recombination coefficients
p.taun = [1e6, 1e6, 1e6];       % [s] SRH time constant for electrons
p.taup = [1e6, 1e6, 1e6];      % [s] SRH time constant for holes

% Raditative recombination could also be set to low values initially if required. 
% p.krad = 1e-20;
% p.kradetl = 1e-20;
% p.kradhtl = 1e-20;

%% General initial parameters
p.tmesh_type = 2;
p.tpoints = 40;

p.Ana = 0;
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
p.mue = [0, 0, 0];          % electron mobility
p.muh = [0, 0, 0];      % hole mobility
p.mui = 0;

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
p.mue = [1e-2, 1e-2, 1e-2];          % electron mobility
p.muh = [1e-2, 1e-2, 1e-2];      % hole mobility
p.mui = 0;

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

p.Ana = 1;
p.calcJ = 0;
p.tmax = 1e-2;
p.t0 = p.tmax/1e10;

soleq.eq = pindrift(sol, p);
disp('Complete')

% 
% %% Set up solution for open circuit
% disp('Switching boundary conditions to zero flux')
% %p.Ana = 0;
% p.BC = 0;
% p.tmax = 1e-9;
% p.t0 = p.tmax/1e3;
% 
% sol = pindrift(soleq.eq, p);
% disp('Complete')
% 
% %% Symmetricise the solution
% disp('Symmetricise solution for open circuit')
% symsol = symmetricize(sol);
% disp('Complete')
% 
% %% Equilibrium solution with mirrored cell and OC boundary conditions, mobility zero
% disp('Initial equilibrium open circuit solution')
% p.BC = BC;
% p.OC = 1;
% p.calcJ = 0;

% %% Switch off mobilities
% p.mue_i = 0;          % electron mobility
% p.muh_i = 0;      % hole mobility
% p.mue_p = 0;
% p.muh_p = 0;
% p.mue_n = 0;
% p.muh_n = 0;
% p.mui = 0;
% 
% ssol = pindrift(symsol, p);
% disp('Complete')
% 
% %% OC with mobility switched on
% disp('Open circuit solution with mobility switched on')
% p.tmax = 1e-6;
% p.t0 = p.tmax/1e3;
% % Switch on mobilities
% p.mue_i = mue_i;          % electron mobility
% p.muh_i = muh_i;      % hole mobility
% p.mue_p = mue_p;
% p.muh_p = muh_p;
% p.mue_n = mue_n;
% p.muh_n = muh_n;
% p.mui = 0;
% 
% ssol = pindrift(ssol, p);
% 
% % Longer time step to ensure equilibrium has been reached
% p.tmax = 1e-2;
% p.t0 = p.tmax/1e3;
% 
% ssol_eq = pindrift(ssol, p);
% disp('Complete')




%% Equilibrium solutions with ion mobility switched on
%% Closed circuit conditions
disp('Closed circuit equilibrium with ions')

p.OC = 0;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p.mui = 1e-6;           % Ions are accelerated to reach equilibrium

sol = pindrift(soleq.eq, p);

% Much longer second step to ensure that ions have migrated
p.calcJ = 0;
p.tmax = 1e-2;
p.t0 = p.tmax/1e3;

soleq.i = pindrift(sol, p);
disp('Complete')

%% Ion equilibrium with surface recombination
disp('Switching on surface recombination')
p.taun = p_original.taun;
p.taup = p_original.taup;

p.calcJ = 0;
p.tmax = 1e-6;
p.t0 = p.tmax/1e3;

soleq.i_sr = pindrift(soleq.i, p);
disp('Complete')

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

disp('EQUILIBRATION COMPLETE')
toc

end
