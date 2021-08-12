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
x_ihalf = par.x_ihalf;
x = xmesh;

%% Time mesh
t = meshgen_t(par);

%% Dependent properties: Prevents recalculation of dependent properties by pdepe defined in Methods
% Can also use AbortSet in class def
Vbi = par.Vbi;
nleft = par.nleft;
nright = par.nright;
pleft = par.pleft;
pright = par.pright;
dev = par.dev;

%% Constants
kB = par.kB;
q = par.q;
e = par.e;
epp0 = par.epp0;

%% Device parameters
device = par.dev_ihalf;
T = par.T;
mue = device.mue;           % Electron mobility
muh = device.muh;           % Hole mobility
Dn = mue*kB*T;              % Electron diffusion coefficient
Dp = muh*kB*T;              % Hole diffusion coefficient
mucat = device.mucat;       % Cation mobility
muani = device.muani;       % Anion mobility
Nc = device.Nc;             % Conduction band effective density of states
Nv = device.Nv;             % Valence band effective density of states
DOScat = device.DOScat;     % Cation density upper limit
DOSani = device.DOSani;     % Anion density upper limit
gradNc = device.gradNc;     % Conduction band effective density of states gradient
gradNv = device.gradNv;     % Valence band effective density of states gradient
gradEA = device.gradEA;     % Electron Affinity gradient
gradIP = device.gradIP;     % Ionisation Potential gradient
epp = device.epp;           % Dielectric constant
eppmax = max(par.epp);      % Maximum dielectric constant (for normalisation)
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
Ncat = device.Ncat;         % Uniform cation density
Nani = device.Nani;         % Uniform anion density 
xprime_n = device.xprime_n;         % Translated x co-ordinates for interfaces
xprime_p = device.xprime_p;         % Translated x co-ordinates for interfaces
sign_xn = device.sign_xn;           % 1 if xn increasing, -1 if decreasing wrt x
sign_xp = device.sign_xp;           % 1 if xp increasing, -1 if decreasing wrt x
alpha0_xn = device.alpha0_xn;       % alpha0_xn is alpha for F = 0 reference to xprime_n
beta0_xp = device.beta0_xp;         % beta0_xp is beta for F = 0 referenced to xprime_p
N_ionic_species = par.N_ionic_species;
N_variables = par.N_ionic_species + 3;  % +3 for V, n, and p
nleft = par.nleft;
nright = par.nright;
pleft = par.pleft;
pright = par.pright;
sn_l = par.sn_l;
sp_l = par.sp_l;
sn_r = par.sn_r;
sp_r = par.sp_r;
Rs = par.Rs;
Rs_initial = par.Rs_initial;

%% Switches and accelerator coefficients
mobset = par.mobset;        % Electronic carrier transport switch
mobseti = par.mobseti;      % Ionic carrier transport switch
K_cation = par.K_cation;    % Cation transport rate multiplier
K_anion = par.K_anion;      % Anion transport rate multiplier
radset = par.radset;        % Radiative recombination switch
SRHset = par.SRHset;        % SRH recombination switch
vsr_zone = device.vsr_zone;
srh_zone = device.vsr_zone;
Field_switch = device.Field_switch;

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

%% Voltage function
Vapp_fun = fun_gen(par.V_fun_type);
Vapp = 0;
Vres = 0;
J = 0;

%% Solver variables
V = 0; n = 0; p = 0; a = 0; c = 0;
dVdx = 0; dndx = 0; dpdx = 0; dadx = 0; dcdx = 0;
C_V = 0; C_n = 0; C_p = 0; C_c = 0; C_a = 0;
F_V = 0; F_n = 0; F_p = 0; F_c = 0; F_a = 0;
S_V = 0; S_n = 0; S_p = 0; S_c = 0; S_a = 0;
r_rad = 0; r_srh = 0; r_vsr = 0; r_np = 0;
alpha = 0; beta = 0;

%% Solver options
% MaxStep = limit maximum time step size during integration
options = odeset('MaxStep', par.MaxStepFactor*0.1*abs(par.tmax - par.t0), 'RelTol', par.RelTol, 'AbsTol', par.AbsTol, 'NonNegative', [1,1,1,0]);

%% Call solver
% inputs with '@' are function handles to the subfunctions
% below for the: equation, initial conditions, boundary conditions
u = pdepe(par.m,@dfpde,@dfic,@dfbc,x,t,options);

%% Ouputs
% Store final voltage reading
par.Vapp = Vapp;

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
        % Get position point
        i = find(x_ihalf <= x);
        i = i(end);

        switch g1_fun_type
            case 'constant'
                gxt1 = int1*gx1(i);
            otherwise
                gxt1 = g1_fun(g1_fun_arg, t)*gx1(i);
        end

        switch g2_fun_type
            case 'constant'
                gxt2 = int2*gx2(i);
            otherwise
                gxt2 = g2_fun(g2_fun_arg, t)*gx2(i);
        end
        g = gxt1 + gxt2;

        %% Unpack Variables       
        switch N_ionic_species
            case 0
                V = u(1);
                n = u(2);
                p = u(3);
                c = 0;
                a = 0;
                dVdx = dudx(1);
                dndx = dudx(2);
                dpdx = dudx(3);
                dcdx = 0;
                dadx = 0;
            case 1
                V = u(1);
                n = u(2);
                p = u(3);
                c = u(4);
                a = Nani(i);
                dVdx = dudx(1);
                dndx = dudx(2);
                dpdx = dudx(3);
                dcdx = dudx(4);
                dadx = 0;
            case 2
                V = u(1);
                n = u(2);
                p = u(3);
                c = u(4);
                a = u(5);
                dVdx = dudx(1);
                dndx = dudx(2);
                dpdx = dudx(3);
                dcdx = dudx(4);
                dadx = dudx(5);
        end
        
        % Volumetric surface recombination gradients
        alpha = sign_xn(i)*q*dVdx/(kB*T) + alpha0_xn(i);
        beta = sign_xp(i)*q*-dVdx/(kB*T) + beta0_xp(i);
        
        %% Equation editor
        % Time-dependence prefactor term
        C_V = 0;
        C_n = 1;
        C_p = 1;
        C_c = 1;
        C_a = 1;
        C = [C_V; C_n; C_p; C_c; C_a];
        
        % Flux terms
        F_V = (epp(i)/eppmax)*dVdx;
        F_n = mue(i)*n*(-dVdx + gradEA(i)) + (Dn(i)*(dndx - ((n/Nc(i))*gradNc(i))));
        F_p = muh(i)*p*(dVdx - gradIP(i)) + (Dp(i)*(dpdx - ((p/Nv(i))*gradNv(i))));
        F_c = mucat(i)*(c*dVdx + kB*T*(dcdx + (c*(dcdx/(DOScat(i)-c)))));
        F_a = muani(i)*(a*-dVdx + kB*T*(dadx+(a*(dadx/(DOSani(i)-a)))));
        F = [F_V; mobset*F_n; mobset*F_p; K_cation*mobseti*F_c; K_anion*mobseti*F_a];
        
        % Electron and holes recombination
        % Radiative
        r_rad = radset*B(i)*((n*p)-(ni(i)^2));
        % Bulk SRH
        r_srh = SRHset*srh_zone(i)*(((n*p)-ni(i)^2)/(taun(i)*(p+pt(i)) + taup(i)*(n+nt(i))));
        % Volumetric surface recombination
        r_vsr = SRHset*vsr_zone(i)*((n*exp(-alpha*xprime_n(i))*p*exp(-beta*xprime_p(i)) - ni(i)^2)...
            /(taun_vsr(i)*(p*exp(-beta*xprime_p(i))+pt(i)) + taup_vsr(i)*(n*exp(-alpha*xprime_n(i))+nt(i))));
        r_np = r_rad + r_srh + r_vsr;
        
        % Source terms
        S_V = Field_switch(i)*(q/(eppmax*epp0))*(-n+p-NA(i)+ND(i)-a+c);
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
    end

%% Define initial conditions.
    function u0 = dfic(x)
        if dficAnalytical
            i = find(xmesh <= x);
            i = i(end);
            
            switch N_ionic_species
                case 0
                    if length(par.dcell) == 1
                        % Single layer
                        u0 = [(x/xmesh(end))*Vbi;
                            nleft*exp((x*(log(nright)-log(nleft)))/par.dcum0(end));
                            pleft*exp((x*(log(pright)-log(pleft)))/par.dcum0(end))];
                    else
                        % Multi-layered
                        u0 = [(x/xmesh(end))*Vbi;
                            dev.n0(i);
                            dev.p0(i)];
                    end

                case 1
                    if length(par.dcell) == 1
                        % Single layer
                        u0 = [(x/xmesh(end))*Vbi;
                            nleft*exp((x*(log(nright)-log(nleft)))/par.dcum0(end));
                            pleft*exp((x*(log(pright)-log(pleft)))/par.dcum0(end));
                            dev.Ncat(i);];
                    else
                        % Multi-layered
                        u0 = [(x/xmesh(end))*Vbi;
                            dev.n0(i);
                            dev.p0(i);
                            dev.Ncat(i);];
                    end
                case 2
                    if length(par.dcell) == 1
                        % Single layer
                        u0 = [(x/xmesh(end))*Vbi;
                            nleft*exp((x*(log(nright)-log(nleft)))/par.dcum0(end));
                            pleft*exp((x*(log(pright)-log(pleft)))/par.dcum0(end));
                            dev.Ncat(i);
                            dev.Nani(i);];
                    else
                        % Multi-layered
                        u0 = [(x/xmesh(end))*Vbi;
                            dev.n0(i);
                            dev.p0(i);
                            dev.Ncat(i);
                            dev.Nani(i);];
                    end
            end
        else
                    u0 = interp1(icx,squeeze(icsol(end,:,:)),x)';
        end
    end

%% Boundary conditions
% Define boundary condtions, refer PDEPE help for the precise meaning of p
% and you l and r refer to left and right.
% in this example I am controlling the flux through the boundaries using
% the difference in concentration from equilibrium and the extraction
% coefficient.
    function [Pl,Ql,Pr,Qr] = dfbc(xl,ul,xr,ur,t)

        switch par.V_fun_type
            case 'constant'
                Vapp = par.V_fun_arg(1);
            otherwise
                Vapp = Vapp_fun(par.V_fun_arg, t);
        end

        switch par.BC
            case 2
                % Non- selective contacts - fixed charge densities for majority carrier
                % and flux for minority carriers- use recombination
                % coefficients sn_l & sp_r to set the surface recombination velocity.
                Pl = [-ul(1);
                    -sn_l*(ul(2) - nleft);
                    ul(3) - pleft;];

                Ql = [0; 1; 0;];

                Pr = [-ur(1)+Vbi-Vapp;
                    ur(2) - nright;
                    sp_r*(ur(3) - pright);];

                Qr = [0; 0; 1;];

                if N_ionic_species == 1 || N_ionic_species == 2
                    % Second element are the boundary condition
                    % coefficients for cations
                    Pl = [Pl; 0];
                    Ql = [Ql; 1];
                    Pr = [Pr; 0];
                    Qr = [Qr; 1];
                end

                if N_ionic_species == 2
                    % Second element are the boundary condition coefficients for anions
                    Pl = [Pl; 0];
                    Ql = [Ql; 1];
                    Pr = [Pr; 0];
                    Qr = [Qr; 1];
                end

            case 3
                % Flux boundary conditions for both carrier types.

                % Calculate series resistance voltage Vres
                if Rs == 0
                    Vres = 0;
                else
                    J = e*sp_r*(ur(3) - pright) - e*sn_r*(ur(2) - nright);
                    if Rs_initial
                        Vres = -J*Rs*t/par.tmax;    % Initial linear sweep
                    else
                        Vres = -J*Rs;
                    end
                end

                Pl = [-ul(1);
                    mobset*(-sn_l*(ul(2) - nleft));
                    mobset*(-sp_l*(ul(3) - pleft));];

                Ql = [0; 1; 1;];

                Pr = [-ur(1)+Vbi-Vapp-Vres;
                    mobset*(sn_r*(ur(2) - nright));
                    mobset*(sp_r*(ur(3) - pright));];

                Qr = [0; 1; 1;];

                if N_ionic_species == 1 || N_ionic_species == 2
                    % Final elements are the boundary condition
                    % coefficients for CATIONS
                    Pl = [Pl; 0];
                    Ql = [Ql; 1];
                    Pr = [Pr; 0];
                    Qr = [Qr; 1];
                end

                if N_ionic_species == 2
                    % Final elements are the boundary condition
                    % coefficients for ANIONS
                    Pl = [Pl; 0];
                    Ql = [Ql; 1];
                    Pr = [Pr; 0];
                    Qr = [Qr; 1];
                end
        end
    end

end
