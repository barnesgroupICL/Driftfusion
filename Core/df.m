function solstruct = df(varargin)
% Core DRIFTFUSION function- organises properties and inputs for pdepe
% A routine to test solving the diffusion and drift equations using the
% matlab pepde solver.
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Solution outputs
% V = u(1) = electrostatic potential
% n = u(2) = electron density
% p = u(3) = holes density
% c = u(4) = cation density (optional)
% a = u(5) = anion density (optional)
%
%% Start code
%% Deal with input arguments
if length(varargin) == 0
    % If no input parameter set then call pc directly
    par = pc;
    dficAnalytical = true;
elseif length(varargin) == 1
    % If one input argument then assume it is the Initial Conditions (IC) solution
    icsol = varargin{1, 1}.u;
    icx = varargin{1, 1}.x;
    par = varargin{1, 1}.par;
    dficAnalytical = false;
elseif length(varargin) == 2
    if max(max(max(varargin{1, 1}.u))) == 0
        par = varargin{2};
        dficAnalytical = true;
    elseif isa(varargin{2}, 'char') == 1            % Checks to see if argument is a character
        input_solstruct = varargin{1, 1};
        par = input_solstruct.par;
        icsol = input_solstruct.u;
        icx = input_solstruct.x;
        dficAnalytical = false;
    else
        input_solstruct = varargin{1, 1};
        icsol = input_solstruct.u;
        icx = input_solstruct.x;
        par = varargin{2};
        dficAnalytical = false;
    end
end

%% Unpack properties
%% Spatial mesh
xmesh = par.xx;
x_sub = par.x_sub;
x = xmesh;

%% Time mesh
t = meshgen_t(par);

%% Dependent properties: Prevents recalculation of dependent properties by pdepe defined in Methods
% Can also use AbortSet in class def
Vbi = par.Vbi;
n0_l = par.n0_l;
n0_r = par.n0_r;
p0_l = par.p0_l;
p0_r = par.p0_r;
dev = par.dev;

%% Constants
kB = par.kB;
q = par.q;
e = par.e;
epp0 = par.epp0;

%% Device parameters
N_ionic_species = par.N_ionic_species;  % Number of ionic species in this solution
N_variables = par.N_ionic_species + 3;  % Number of variables in this solution (+3 for V, n, and p)
N_max_variables = par.N_max_variables;  % Maximum number of variables in this version

device = par.dev_sub;
T = par.T;
mu_n = device.mu_n;         % Electron mobility
mu_p = device.mu_p;         % Hole mobility
mu_c = device.mu_c;         % Cation mobility
mu_a = device.mu_a;         % Anion mobility
Nc = device.Nc;             % Conduction band effective density of states
Nv = device.Nv;             % Valence band effective density of states
c_max = device.c_max;       % Cation density upper limit
a_max = device.a_max;       % Anion density upper limit
gradNc = device.gradNc;     % Conduction band effective density of states gradient
gradNv = device.gradNv;     % Valence band effective density of states gradient
gradEA = device.gradEA;     % Electron Affinity gradient
gradIP = device.gradIP;     % Ionisation Potential gradient
epp = device.epp;           % Dielectric constant
epp_factor = par.epp_factor;      % Maximum dielectric constant (for normalisation)
B = device.B;               % Radiative recombination rate coefficient
ni = device.ni;             % Intrinsic carrier density
taun = device.taun;         % Electron SRH time constant
taup = device.taup;         % Electron SRH time constant
taun_vsr = device.taun_vsr; % Electron SRH time constant- volumetric interfacial surface recombination scheme
taup_vsr = device.taup_vsr; % Electron SRH time constant- volumetric interfacial surface recombination scheme
nt = device.nt;             % SRH electron trap constant
pt = device.pt;             % SRH hole trap constant
NA = device.NA;             % Acceptor doping density
ND = device.ND;             % Donor doping density
% Set up counter ion density arrays
switch N_ionic_species
    case 0                  % Nani, Ncat, a, and c set to zero for Poisson
        Ncat = zeros(1, length(x_sub));
        Nani = zeros(1, length(x_sub));
    case 1                  % Nani and a both set to zero for Poisson
        Ncat = device.Ncat;
        Nani = zeros(1, length(x_sub));
    case 2
        Ncat = device.Ncat;
        Nani = device.Nani;
end
xprime_n = device.xprime_n;         % Translated x co-ordinates for interfaces
xprime_p = device.xprime_p;         % Translated x co-ordinates for interfaces
sign_xn = device.sign_xn;           % 1 if xn increasing, -1 if decreasing wrt x
sign_xp = device.sign_xp;           % 1 if xp increasing, -1 if decreasing wrt x
alpha0_xn = device.alpha0_xn;       % alpha0_xn is alpha for F = 0 reference to xprime_n
beta0_xp = device.beta0_xp;         % beta0_xp is beta for F = 0 referenced to xprime_p

z_c = par.z_c;
z_a = par.z_a;
n0_l = par.n0_l;
n0_r = par.n0_r;
p0_l = par.p0_l;
p0_r = par.p0_r;
sn_l = par.sn_l;
sp_l = par.sp_l;
sn_r = par.sn_r;
sp_r = par.sp_r;
Rs = par.Rs;
gamma = par.gamma;          % Blakemore approximation coefficient, 0 for Boltzmann stats

%% Switches and accelerator coefficients
mobset = par.mobset;        % Electronic carrier transport switch
mobseti = par.mobseti;      % Ionic carrier transport switch
K_c = par.K_c;              % Cation transport rate multiplier
K_a = par.K_a;              % Anion transport rate multiplier
radset = par.radset;        % Radiative recombination switch
SRHset = par.SRHset;        % SRH recombination switch
vsr_zone = device.vsr_zone;
srh_zone = device.srh_zone;
Rs_initial = par.Rs_initial;
Field_switch = dev.Field_switch;

%% Generation
g1_fun = fun_gen(par.g1_fun_type);
g2_fun = fun_gen(par.g2_fun_type);
gxt1 = 0;
gxt2 = 0;
g = 0;
gx1 = par.gx1;
gx2 = par.gx2;
int1 = par.int1;
int2 = par.int2;
g1_fun_type = par.g1_fun_type;
g2_fun_type = par.g2_fun_type;
g1_fun_arg = par.g1_fun_arg;
g2_fun_arg = par.g2_fun_arg;

% Illumination type g1_fun_type and g2_fun_type - convert to Boolean for
% faster execution in PDEPE
if strcmp(g1_fun_type, 'constant')
    g1_fun_type_constant = 1;
else
    g1_fun_type_constant = 0;
end
if strcmp(g2_fun_type, 'constant')
    g2_fun_type_constant = 1;
else
    g2_fun_type_constant = 0;
end

% Check for negative generation and deal error if present
gM = g1_fun(g1_fun_arg, t')*gx1 + g2_fun(g2_fun_arg, t')*gx2;
if any(any(gM < 0))
    error('Generation cannot be negative - please check your generation function and associated inputs')
end

%% Voltage function
Vapp_fun = fun_gen(par.V_fun_type);
Vres = 0;
J = 0;

%% Solver variables
i = 1;
V = 0; n = 0; p = 0; a = 0; c = 0;
dVdx = 0; dndx = 0; dpdx = 0; dadx = 0; dcdx = 0;
F_V = 0; F_n = 0; F_p = 0; F_c = 0; F_a = 0;
S_V = 0; S_n = 0; S_p = 0; S_c = 0; S_a = 0;
r_rad = 0; r_srh = 0; r_vsr = 0; r_np = 0;
alpha = 0; beta = 0;
G_n = 1;    % Diffusion enhancement prefactor electrons
G_p = 1;    % Diffusion enhancement prefactor holes

%% Initialise solution arrays
u_maxvar = zeros(N_max_variables, 1);
dudx_maxvar = zeros(N_max_variables, 1);
ul_maxvar = zeros(N_max_variables, 1);
ur_maxvar = zeros(N_max_variables, 1);

%% Solver options
% MaxStep = limit maximum time step size during integration
options = odeset('MaxStep', par.MaxStepFactor*0.1*par.tmax, 'RelTol', par.RelTol, 'AbsTol', par.AbsTol);

%% Call solver
% inputs with '@' are function handles to the subfunctions
% below for the: equation, initial conditions, boundary conditions
u = pdepe(par.m,@dfpde,@dfic,@dfbc,x,t,options);

%% Outputs
% Solutions and meshes to structure
solstruct.u = u;
solstruct.x = x;
solstruct.t = t;

% Store parameters object
solstruct.par = par;

%% Volumetric surface recombination error check
if par.vsr_mode == 1 && par.vsr_check == 1
    compare_rec_flux(solstruct, par.RelTol_vsr, par.AbsTol_vsr, 0);
end

%% Subfunctions
% Set up partial differential equation (pdepe) (see MATLAB pdepe help for details of C,F,S)
% C = Time-dependence prefactor
% F = Flux terms
% S = Source terms
    function [C,F,S] = dfpde(x,t,u,dudx)
        % reset position point
        if x == x_sub(1)
            i = 1;
        end
        
        % Generation
        if g1_fun_type_constant
            gxt1 = int1*gx1(i);
        else
            gxt1 = g1_fun(g1_fun_arg, t)*gx1(i);
        end
        
        if g2_fun_type_constant
            gxt2 = int2*gx2(i);
        else
            gxt2 = g2_fun(g2_fun_arg, t)*gx2(i);
        end
        g = gxt1 + gxt2;
        
        %% Unpack Variables
        u_maxvar(1:N_variables) = u;
        dudx_maxvar(1:N_variables) = dudx;
        V = u_maxvar(1);
        n = u_maxvar(2);
        p = u_maxvar(3);
        c = u_maxvar(4);
        a = u_maxvar(5);
        dVdx = dudx_maxvar(1);
        dndx = dudx_maxvar(2);
        dpdx = dudx_maxvar(3);
        dcdx = dudx_maxvar(4);
        dadx = dudx_maxvar(5);
        
        % Diffusion enhancement prefactors (gamma = 0 for Boltz)
        G_n = Nc(i)/(Nc(i) - gamma*n);
        G_p = Nv(i)/(Nv(i) - gamma*p);
        
        %% Equation editor
        % Time-dependence pre-factor (pre-allocated above)
        % Time-dependence prefactor term
        C_V = 0;
        C_n = 1;
        C_p = 1;
        C_c = 1;
        C_a = 1;
        C = [C_V; C_n; C_p; C_c; C_a];
        
        % Flux terms
        F_V = (epp(i)/epp_factor)*dVdx;
        F_n = mu_n(i)*n*(-dVdx + gradEA(i)) + (G_n*mu_n(i)*kB*T*(dndx - ((n/Nc(i))*gradNc(i))));
        F_p = mu_p(i)*p*(dVdx - gradIP(i)) + (G_p*mu_p(i)*kB*T*(dpdx - ((p/Nv(i))*gradNv(i))));
        F_c = mu_c(i)*(z_c*c*dVdx + kB*T*(dcdx + (c*(dcdx/(c_max(i) - c)))));
        F_a = mu_a(i)*(z_a*a*dVdx + kB*T*(dadx + (a*(dadx/(a_max(i) - a)))));
        F = [F_V; mobset*F_n; mobset*F_p; mobseti*K_c*F_c; mobseti*K_a*F_a];
        
        % Electron and hole recombination
        % Radiative
        r_rad = radset*B(i)*(n*p - ni(i)^2);
        % Bulk SRH
        r_srh = SRHset*srh_zone(i)*((n*p - ni(i)^2)/(taun(i)*(p + pt(i)) + taup(i)*(n + nt(i))));
        % Volumetric surface recombination
        alpha = (sign_xn(i)*q*dVdx/(kB*T)) + alpha0_xn(i);
        beta = (sign_xp(i)*q*-dVdx/(kB*T)) + beta0_xp(i);
        r_vsr = SRHset*vsr_zone(i)*((n*exp(-alpha*xprime_n(i))*p*exp(-beta*xprime_p(i)) - ni(i)^2)...
            /(taun_vsr(i)*(p*exp(-beta*xprime_p(i)) + pt(i)) + taup_vsr(i)*(n*exp(-alpha*xprime_n(i)) + nt(i))));
        % Total electron and hole recombination
        r_np = r_rad + r_srh + r_vsr;
        
        % Source terms
        S_V = (1/(epp_factor*epp0))*(-n + p - NA(i) + ND(i) + z_a*a + z_c*c - z_a*Nani(i) - z_c*Ncat(i));
        S_n = g - r_np;
        S_p = g - r_np;
        S_c = 0;
        S_a = 0;
        S = [S_V; S_n; S_p; S_c; S_a];
        
        % Remove unused variables - faster and tidier than using conditional
        % statements
        C = C(1:N_variables);
        F = F(1:N_variables);
        S = S(1:N_variables);
        
        i = i+1;
    end

%% Initial conditions.
    function u0 = dfic(x)
        if x == x_sub(1)
            i = 1;
        end
        
        if length(par.dcell) == 1
            % Single layer
            u0_ana = [(x/xmesh(end))*Vbi;
                n0_l*exp((x*(log(n0_r)-log(n0_l)))/par.dcum0(end));
                p0_l*exp((x*(log(p0_r)-log(p0_l)))/par.dcum0(end));
                dev.Ncat(i);
                dev.Nani(i);];
        else
            % Multi-layered
            u0_ana = [(x/xmesh(end))*Vbi;
                dev.n0(i);
                dev.p0(i);
                dev.Ncat(i);
                dev.Nani(i);];
        end
        u0_ana = u0_ana(1:N_variables);
        
        % Organise ICs based on number of variables and SOL_IC
        if dficAnalytical
            u0 = u0_ana;
        else
            % Initial conditions taken from input solution
            u0_input = interp1(icx,squeeze(icsol(end,:,:)),x)';
            % If the number of variables has increased then add analytical
            % ICs for missing variables
            if N_variables > length(u0_input)
                % Add initial conditions for new variables from U_ANA
                u0(1:length(u0_input), 1) = u0_input;
                u0(length(u0_input)+1:N_variables, 1) = u0_ana(length(u0_input)+1:N_variables);
            else
                u0 = u0_input;
            end
        end
        
        i = i+1;
    end

%% Boundary conditions
% Refer to PDEPE help for the precise meaning of P and Q
% l and r refer to left and right boundaries.
    function [Pl,Ql,Pr,Qr] = dfbc(xl,ul,xr,ur,t)
        
        ul_maxvar(1:N_variables) = ul;
        ur_maxvar(1:N_variables) = ur;
        
        V_l = ul_maxvar(1);
        V_r = ur_maxvar(1);
        n_l = ul_maxvar(2);
        n_r = ur_maxvar(2);
        p_l = ul_maxvar(3);
        p_r = ur_maxvar(3);
        c_l = ul_maxvar(4);
        c_r = ur_maxvar(4);
        a_l = ul_maxvar(5);
        a_r = ur_maxvar(5);
        
        switch par.V_fun_type
            case 'constant'
                Vapp = par.V_fun_arg(1);
            otherwise
                Vapp = Vapp_fun(par.V_fun_arg, t);
        end
        
        % Flux boundary conditions for both carrier types.
        % Calculate series resistance voltage Vres
        if Rs == 0
            Vres = 0;
        else
            J = e*sp_r*(p_r - p0_r) - e*sn_r*(n_r - n0_r);
            if Rs_initial
                Vres = -J*Rs*t/par.tmax;    % Initial linear sweep
            else
                Vres = -J*Rs;
            end
        end
        
        Pl = [-V_l;
            mobset*(-sn_l*(n_l - n0_l));
            mobset*(-sp_l*(p_l - p0_l));
            0;
            0;];
        
        Ql = [0;
            1;
            1;
            1;
            1;];
        
        Pr = [-V_r+Vbi-Vapp-Vres;
            mobset*(sn_r*(n_r - n0_r));
            mobset*(sp_r*(p_r - p0_r));
            0;
            0;];
        
        Qr = [0;
            1;
            1;
            1;
            1;];
        
        % Remove unused entries
        Pl = Pl(1:N_variables);
        Pr = Pr(1:N_variables);
        Ql = Ql(1:N_variables);
        Qr = Qr(1:N_variables);
    end
end
