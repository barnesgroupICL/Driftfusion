function [sol_eq, sol_i_eq, ssol_eq, ssol_i_eq, sol_i_eq_SR, ssol_i_eq_SR, ssol_i_1S_SR ] = equilibrate(option)
% Uses analytical initial conditions and runs to equilibrium
% Note that tmax is consistently adjusted to appropriate values for to
% ensure there are numerous mesh points where large gradients in the time
% dimension are present

%% options
% 1 = closed circuit only no interfacial rec
% 2 = closed circuit including interfacial rec
% 3 = open circuit

tic;    % Start stopwatch

%% Initial arguments
% Setting sol.sol = 0 enables a parameters structure to be read into
% pindrift but indicates that the initial conditions should be the
% analytical solutions
sol.sol = 0;    

p = pinParams;
p_initial = p;       % Store initial params

%% Start with low recombination coefficients
p.klin = 0;
p.klincon = 0;
p.taun_etl = 1e6;       % [s] SRH time constant for electrons
p.taup_etl = 1e6;      % [s] SRH time constant for holes
p.taun_htl = 1e6;       %%%% USE a high value of (e.g.) 1 to switch off
p.taup_htl = 1e6;

% Raditative recombination could also be set to low values initially if required. 
% p.krad = 1e-20;
% p.kradetl = 1e-20;
% p.kradhtl = 1e-20;

%% General initial parameters
p.tmesh_type = 2;
p.tpoints = 200;

p.Ana = 0;
p.JV = 0;
p.Vapp = 0;
p.Int = 0;
p.pulseon = 0; 
p.OC = 0;
p.p_initial.BC = p_initial.BC;
p.tmesh_type = 2;
p.tmax = 1e-9;
p.t0 = p.tmax/1e4;

%% Switch off mobilities
p.mue_i = 0;          % electron mobility
p.muh_i = 0;      % hole mobility
p.mue_p = 0;
p.muh_p = 0;
p.mue_n = 0;
p.muh_n = 0;
p.mui = 0;

% Switch off extraction and recombination
p.sn_ext = 0;
p.sn_rec = 0;
p.sp_ext = 0;
p.sp_rec = 0;

%% Initial solution with zero mobility
disp('Initial solution, zero mobility')
sol = pindrift(sol, p);
disp('Complete')

% Switch on mobilities
p.mue_i = p_initial.mue_i;          % electron mobility
p.muh_i = p_initial.muh_i;      % hole mobility
p.mue_p = p_initial.mue_p;
p.muh_p = p_initial.muh_p;
p.mue_n = p_initial.mue_n;
p.muh_n = p_initial.muh_n;
p.mui = 0;

p.sn_ext = p_initial.sn_ext;
p.sn_rec = p_initial.sn_rec;
p.sp_ext = p_initial.sp_ext;
p.sp_rec = p_initial.sp_rec;

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

sol_eq = pindrift(sol, p);
disp('Complete')

if option == 3
    
    %% Set up solution for open circuit
    disp('Switching boundary conditions to zero flux')
    %p.Ana = 0;
    p.p_initial.BC = 0;
    p.tmax = 1e-9;
    p.t0 = p.tmax/1e3;
    
    sol = pindrift(sol_eq, p);
    disp('Complete')
    
    %% Symmetricise the solution
    disp('Symmetricise solution for open circuit')
    symsol = symmetricize(sol);
    disp('Complete')

    %% Equilibrium solution with mirrored cell and OC boundary conditions, mobility zero
    disp('Initial equilibrium open circuit solution')
    p.p_initial.BC = p_initial.BC;
    p.OC = 1;
    p.calcJ = 0;

    %% Switch off mobilities
    p.mue_i = 0;          % electron mobility
    p.muh_i = 0;      % hole mobility
    p.mue_p = 0;
    p.muh_p = 0;
    p.mue_n = 0;
    p.muh_n = 0;
    p.mui = 0;

    ssol = pindrift(symsol, p);
    disp('Complete')

    %% OC with mobility switched on
    disp('Open circuit solution with mobility switched on')
    p.tmax = 1e-6;
    p.t0 = p.tmax/1e3;
    % Switch on mobilities
    p.mue_i = p_initial.mue_i;          % electron mobility
    p.muh_i = p_initial.muh_i;      % hole mobility
    p.mue_p = p_initial.mue_p;
    p.muh_p = p_initial.muh_p;
    p.mue_n = p_initial.mue_n;
    p.muh_n = p_initial.muh_n;
    p.mui = 0;

    ssol = pindrift(ssol, p);

    % Longer time step to ensure equilibrium has been reached
    p.tmax = 1e-2;
    p.t0 = p.tmax/1e3;

    ssol_eq = pindrift(ssol, p);
    disp('Complete')

end

%% Equilibrium solutions with ion mobility switched on
%% Closed circuit conditions
disp('Closed circuit equilibrium with ions')

p.OC = 0;
p.tmax = 1e-9;
p.t0 = p.tmax/1e3;
p.mui = 1e-6;           % Ions are accelerated to reach equilibrium

sol = pindrift(sol_eq, p);

% Much longer second step to ensure that ions have migrated
p.calcJ = 2;
p.tmax = 1e-2;
p.t0 = p.tmax/1e3;

sol_i_eq = pindrift(sol, p);
disp('Complete')

if option == 2 || option == 3
    
    %% Ion equilibrium with surface recombination
    disp('Switching on surface recombination')
    p.taun_etl = p_initial.taun_etl;
    p.taup_etl = p_initial.taup_etl;
    p.taun_htl = p_initial.taun_htl;
    p.taup_htl = p_initial.taup_htl;
    
    p.calcJ = 0;
    p.tmax = 1e-6;
    p.t0 = p.tmax/1e3;
    
    sol_i_eq_SR = pindrift(sol_i_eq, p);
    disp('Complete')
    
    % Switch off SR
    p.taun_etl = 1e6;
    p.taup_etl = 1e6;
    p.taun_htl = 1e6;
    p.taup_htl = 1e6;
    
end

if option == 3
    
    %% Symmetricise closed circuit condition
    disp('Symmetricise equilibriumion solution')
    symsol = symmetricize(sol_i_eq);
    disp('Complete')
    
    p.OC = 1;
    p.tmax = 1e-9;
    p.t0 = p.tmax/1e3;
    
    %% Switch off mbilities
    p.mue_i = 0;          % electron mobility
    p.muh_i = 0;      % hole mobility
    p.mue_p = 0;
    p.muh_p = 0;
    p.mue_n = 0;
    p.muh_n = 0;
    p.mui = 0;
    
    %% OC condition with ions at equilbirium
    disp('Open circuit equilibrium with ions')
    ssol_sym = pindrift(symsol, p);
    
    p.tmax = 1e-9;
    p.t0 = p.tmax/1e3;
    
    % Switch on mobilities
    p.mue_i = p_initial.mue_i;          % electron mobility
    p.muh_i = p_initial.muh_i;      % hole mobility
    p.mue_p = p_initial.mue_p;
    p.muh_p = p_initial.muh_p;
    p.mue_n = p_initial.mue_n;
    p.muh_n = p_initial.muh_n;
    p.mui = 0;
    
    ssol = pindrift(ssol_sym, p);
    
    % Switch on ion mobility to ensure equilibrium has been reached
    p.tmax = 1e-9;
    p.t0 = p.tmax/1e3;
    p.mui = 1e-6;
    
    ssol = pindrift(ssol, p);
    
    p.tmax = 1;
    p.t0 = p.tmax/1e6;
    
    ssol_i_eq = pindrift(ssol, p);
    
    disp('Complete')
    
    %% Ions, OC Surface recombination
    p.taun_etl = p_initial.taun_etl;
    p.taup_etl = p_initial.taup_etl;
    p.taun_htl = p_initial.taun_htl;
    p.taup_htl = p_initial.taup_htl;
    
    p.tmax = 1e-3;
    p.t0 = p.tmax/1e6;
    
    ssol_i_eq_SR = pindrift(ssol_i_eq , p);
    
    %% 1 Sun quasi equilibrium
    disp('1 Sun quasi equilibrium')
    p.tmax = 1e-3;
    p.t0 = p.tmax/1e6;
    p.Int = 1;
    
    ssol_i_1S_SR = pindrift(ssol_i_eq_SR, p);
    
    disp('Complete')
    
end

if option ~= 2
    
    sol_i_eq_SR = 0;
    
end

if option ~= 3
    
    ssol_eq = 0;
    ssol_i_eq = 0;
    ssol_i_eq_SR = 0;
    ssol_i_1S_SR = 0;

disp('EQUILIBRATION COMPLETE')
toc

end
