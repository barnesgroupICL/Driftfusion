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
%% Start code
% n = u(1) = electron density
% p = u(2) = holes density
% c = u(3) = cation density
% V = u(4) = electrostatic potential
% a = u(5) = anion density

%% Deal with input arguments
if length(varargin) == 0
    % If no input parameter set then call pc directly
    par = pc;
elseif length(varargin) == 1
    % If one input argument then assume it is the Initial Conditions (IC) solution
    icsol = varargin{1, 1}.u;
    icx = varargin{1, 1}.x;
    par = pc;
elseif length(varargin) == 2
    if max(max(max(varargin{1, 1}.u))) == 0
        par = varargin{2};
    elseif isa(varargin{2}, 'char') == 1            % Checks to see if argument is a character

        input_solstruct = varargin{1, 1};
        par = input_solstruct.par;
        icsol = input_solstruct.u;
        icx = input_solstruct.x;
    else
        input_solstruct = varargin{1, 1};
        icsol = input_solstruct.u;
        icx = input_solstruct.x;
        par = varargin{2};
    end
end

%% Unpack properties
% Dependent properties: Prevents recalculation of dependent properties by pdepe defined in Methods
% Can also use AbortSet in class def
Vbi = par.Vbi;
nleft = par.nleft;
nright = par.nright;
pleft = par.pleft;
pright = par.pright;
xmesh = par.xx;
x_ihalf = par.x_ihalf;
devihalf = par.dev_ihalf;
dev = par.dev;

%% Constants
kB = par.kB;
T = par.T;
q = par.q;
epp0 = par.epp0;

%% Switches and accelerator coefficients
mobset = par.mobset;        % Electronic carrier transport switch
mobseti = par.mobseti;      % Ionic carrier transport switch
K_cation = par.K_cation;    % Cation transport rate multiplier
K_anion = par.K_anion;      % Anion transport rate multiplier
radset = par.radset;        % Radiative recombination switch
SRHset = par.SRHset;        % SRH recombination switch

%% Device parameters
device = par.dev_ihalf;
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
B = device.B;         % Radiative recombination rate coefficient
ni = device.ni;             % Intrinsic carrier density
n0 = device.n0;             % Intrinsic carrier density
p0 = device.p0;             % Intrinsic carrier density
taun = device.taun;         % Electron SRH time constant
taup = device.taup;         % Electron SRH time constant
nt = device.nt;             % SRH electron trap constant
pt = device.pt;             % SRH hole trap constant
NA = device.NA;             % Acceptor doping density
ND = device.ND;             % Donor doping density
Ncat = device.Ncat;         % Uniform cation density
Nani = device.Nani;         % Uniform anion density
N_ionic_species = par.N_ionic_species;      % Number of ionic species
nleft = par.nleft;
nright = par.nright;
pleft = par.pleft;
pright = par.pright;
sn_l = par.sn_l;
sp_l = par.sp_l;
sn_r = par.sn_r;
sp_r = par.sp_r;

%% Spatial mesh
x = xmesh;

%% Time mesh
t = meshgen_t(par);
tmesh = t;

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
Jr = 0;

%% Solver options
% MaxStep = limit maximum time step size during integration
options = odeset('MaxStep', par.MaxStepFactor*0.1*abs(par.tmax - par.t0), 'RelTol', par.RelTol, 'AbsTol', par.AbsTol, 'NonNegative', [1,1,1,0]);

%% Call solver
% inputs with '@' are function handles to the subfunctions
% below for the: equation, initial conditions, boundary conditions
u = pdepe(par.m,@dfpde,@dfic,@dfbc,x,t,options);

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

        %% Variables
        n = u(1); p = u(2); c = u(3); V = u(4);

        if N_ionic_species == 2
            a = u(5);           % Include anion variable
            dadx = dudx(5);
        else
            a = Nani(i);        % Otherwise set anion density to be background
        end

        %% Gradients
        dndx = dudx(1); dpdx = dudx(2); dcdx = dudx(3); dVdx = dudx(4);

        %% Equation editor
        % Time-dependence prefactor term
        C_electron = 1;
        C_hole = 1;
        C_cation = 1;
        C_potential = 0;
        C = [C_electron; C_hole; C_cation; C_potential];

        % Flux terms
        F_electron   = mue(i)*n*(-dVdx + gradEA(i)) + (Dn(i)*(dndx - ((n/Nc(i))*gradNc(i))));
        F_hole       = muh(i)*p*(dVdx - gradIP(i)) + (Dp(i)*(dpdx - ((p/Nv(i))*gradNv(i))));
        F_cation     = mucat(i)*(c*dVdx + kB*T*(dcdx + (c*(dcdx/(DOScat(i)-c)))));
        F_potential  = (epp(i)/eppmax)*dVdx;
        F = [mobset*F_electron; mobset*F_hole; K_cation*mobseti*F_cation; F_potential];

        % Source terms - simplified first order recombination terms
        S_electron = g - SRHset*(((n-n0(i))/taun(i)) + ((p-p0(i))/taup(i))); %- SRHset*(((n*p)-ni(i)^2)/((taun(i)*(p+pt(i)))+(taup(i)*(n+nt(i)))));
        S_hole     = g - SRHset*(((n-n0(i))/taun(i)) + ((p-p0(i))/taup(i)));%- SRHset*(((n*p)-ni(i)^2)/((taun(i)*(p+pt(i)))+(taup(i)*(n+nt(i)))));
        S_cation   = 0;
        S_potential = (q/(eppmax*epp0))*(-n+p-NA(i)+ND(i)-a+c+Nani(i)-Ncat(i));
        S = [S_electron; S_hole; S_cation; S_potential];

        if N_ionic_species == 2     % Condition for anion terms
            C_anion = 1;
            C = [C; C_anion];

            F_anion = muani(i)*(a*-dVdx + kB*T*(dadx+(a*(dadx/(DOSani(i)-a)))));
            F = [F; K_anion*mobseti*(F_anion)];

            S_anion = 0;
            S = [S; S_anion];
        end
    end

%% Define initial conditions.
    function u0 = dfic(x)

        if isempty(varargin) || length(varargin) >= 1 && max(max(max(varargin{1, 1}.u))) == 0

            i = find(xmesh <= x);
            i = i(end);

            switch N_ionic_species
                case 1
                    if length(par.dcell) == 1
                        % Single layer
                        u0 = [nleft*exp((x*(log(nright)-log(nleft)))/par.dcum0(end));
                            pleft*exp((x*(log(pright)-log(pleft)))/par.dcum0(end));
                            dev.Ncat(i);
                            (x/xmesh(end))*Vbi;];
                    else
                        % Multi-layered
                        u0 = [dev.n0(i);
                            dev.p0(i);
                            dev.Ncat(i);
                            (x/xmesh(end))*Vbi;];
                    end
                case 2
                    if length(par.dcell) == 1
                        % Single layer
                        u0 = [nleft*exp((x*(log(nright)-log(nleft)))/par.dcum0(end));
                            pleft*exp((x*(log(pright)-log(pleft)))/par.dcum0(end));
                            dev.Ncat(i);
                            (x/xmesh(end))*Vbi;
                            dev.Nani(i);];
                    else
                        % Multi-layered
                        u0 = [dev.n0(i);
                            dev.p0(i);
                            dev.Ncat(i);
                            (x/xmesh(end))*Vbi;
                            dev.Nani(i);];
                    end
            end
        elseif length(varargin) == 1 || length(varargin) >= 1 && max(max(max(varargin{1, 1}.u))) ~= 0
            switch par.N_ionic_species
                case 1
                    u0 = [interp1(icx,icsol(end,:,1),x);
                        interp1(icx,icsol(end,:,2),x);
                        interp1(icx,icsol(end,:,3),x);
                        interp1(icx,icsol(end,:,4),x);];
                case 2
                    % insert previous solution and interpolate the x points
                    u0 = [interp1(icx,icsol(end,:,1),x);
                        interp1(icx,icsol(end,:,2),x);
                        interp1(icx,icsol(end,:,3),x);
                        interp1(icx,icsol(end,:,4),x);
                        interp1(icx,icsol(end,:,5),x);];
            end
        end
    end

%% Boundary conditions
% Define boundary condtions, refer PDEPE help for the precise meaning of p
% and you l and r refer to left and right.
% in this example I am controlling the flux through the boundaries using
% the difference in concentration from equilibrium and the extraction
% coefficient.
    function [pl,ql,pr,qr] = dfbc(xl,ul,xr,ur,t)

        switch par.V_fun_type
            case 'constant'
                Vapp = par.V_fun_arg(1);
            otherwise
                Vapp = Vapp_fun(par.V_fun_arg, t);
        end

        % Easy variable names
        nl = ul(1); pl = ul(2); cl = ul(3); Vl = ul(4);
        nr = ur(1); pr = ur(2); cr = ur(3); Vr = ur(4);

        if N_ionic_species == 2
            al = ul(5); ar = ur(5);
        end

        switch par.BC
            case 2
                % Non- selective contacts - fixed charge densities for majority carrier
                % and flux for minority carriers- use recombination
                % coefficients sn_l & sp_r to set the surface recombination velocity.
                pl = [-sn_l*(ul(1) - nleft);
                    ul(2) - pleft;
                    0;
                    -ul(4);];

                ql = [1; 0;1; 0;];

                pr = [ur(1) - nright;
                    sp_r*(ur(2) - pright);
                    0;
                    -ur(4)+Vbi-Vapp;];

                qr = [0; 1; 1; 0;];

                if N_ionic_species == 2
                    % Second element are the boundary condition coefficients for anions
                    pl = [pl; 0];
                    ql = [ql; 1];
                    pr = [pr; 0];
                    qr = [qr; 1];
                end

            case 3
                % Flux boundary conditions for both carrier types.

                % Calculate series resistance voltage Vres
                if par.Rs == 0
                    Vres = 0;
                else
                    Jr = par.e*sp_r*(ur(2) - pright) - par.e*sn_r*(ur(1) - nright);
                    if par.Rs_initial
                        Vres = Jr*par.Rs*t/par.tmax;    % Initial linear sweep
                    else
                        Vres = Jr*par.Rs;
                    end
                end

                pl = [mobset*(-sn_l*(ul(1) - nleft));
                    mobset*(-sp_l*(ul(2) - pleft));
                    0;
                    -ul(4);];

                ql = [1; 1; 1; 0;];

                pr = [mobset*(sn_r*(ur(1) - nright));
                    mobset*(sp_r*(ur(2) - pright));
                    0;
                    -ur(4)+Vbi-Vapp+Vres;];

                qr = [1; 1; 1; 0;];

                if N_ionic_species == 2
                    %% Second element are the boundary condition coefficients for anions
                    pl = [pl; 0];
                    ql = [ql; 1];
                    pr = [pr; 0];
                    qr = [qr; 1];
                end

        end
    end

%% Ouputs
% Store final voltage reading
par.Vapp = Vapp;

% Solutions and meshes to structure
solstruct.u = u;
solstruct.x = x;
solstruct.t = t;

% Store parameters object
solstruct.par = par;

end
