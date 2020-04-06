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
% V = u(1) = electrostatic potential
% n = u(2) = electron density
% p = u(3) = holes density
% c = u(4) = cation density (optional)
% a = u(5) = anion density (optional)

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

%% Prepare variables for solver

% prepare search for x in x_ihalf
x_ihalf_length = length(x_ihalf);
sampling_step = round(sqrt(x_ihalf_length));
x_ihalf_padded = [x_ihalf, nan(1,sampling_step-mod(x_ihalf_length,sampling_step))];
x_ihalf_padded_mat = reshape(x_ihalf_padded,sampling_step,[]);
x_ihalf_sampled = x_ihalf_padded_mat(1,:);

% normalise epp
epp_norm = epp/eppmax;

% prepare constant generation profiles
int1gx1 = int1*gx1;
int2gx2 = int2*gx2;

% S_potential prefactor
q_over_eppmax_epp0 = q/(eppmax*epp0);

% pre-calculation S_potential static charge
NANDNaniNcat = -NA+ND+Nani-Ncat;

% pre-calculate kB*T
kBT = kB*T;

% convert to boolean for faster check in dfode
g1_fun_type_constant = g1_fun_type == "constant";
g2_fun_type_constant = g2_fun_type == "constant";

% C: Time-dependence prefactor term
switch N_ionic_species
    case 1
        C_potential = 0;
        C_electron = 1;
        C_hole = 1;
        C_cation = 1;
        Cpre = [C_potential; C_electron; C_hole; C_cation];
    case 2
        C_potential = 0;
        C_electron = 1;
        C_hole = 1;
        C_cation = 1;
        C_anion = 1;
        Cpre = [C_potential; C_electron; C_hole; C_cation; C_anion];
    otherwise
        C_potential = 0;
        C_electron = 1;
        C_hole = 1;
        Cpre = [C_potential; C_electron; C_hole];
end

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
    
    % C: Time-dependence prefactor term
    C = Cpre;

    % Get position point
    sampled_i = find(x_ihalf_sampled <= x);
    sampled_i = sampled_i(end);
    i = (sampled_i-1)*sampling_step + find(x == x_ihalf_padded_mat(:,sampled_i));

    % g: Generation terms
    if g1_fun_type_constant
        gxt1 = int1gx1(i);
    else
        gxt1 = g1_fun(g1_fun_arg, t)*gx1(i);
    end

    if g2_fun_type_constant
        gxt2 = int2gx2(i);
    else
        gxt2 = g2_fun(g2_fun_arg, t)*gx2(i);
    end
        
        g = gxt1 + gxt2;
        
        %% Variables
        V = u(1); n = u(2); p = u(3); 

        if N_ionic_species == 1
            c = u(4);           % Include cation variable
            dcdx = dudx(4);
            a = Nani(i);
        elseif N_ionic_species == 2
            c = u(4);           % Include cation variable
            a = u(5);           % Include anion variable
            dcdx = dudx(4);
            dadx = dudx(5);    
        else
            c = Ncat(i);
            a = Nani(i);
        end
        
        %% Gradients
        dVdx = dudx(1); dndx = dudx(2); dpdx = dudx(3); 
        
        % Flux terms
        F_potential  = epp_norm(i)*dVdx;
        F_electron   = mue(i)*n*(-dVdx + gradEA(i)) + (Dn(i)*(dndx - ((n/Nc(i))*gradNc(i))));
        F_hole       = muh(i)*p*(dVdx - gradIP(i)) + (Dp(i)*(dpdx - ((p/Nv(i))*gradNv(i))));      
        F = [F_potential; mobset*F_electron; mobset*F_hole]; 
            
        % Source terms
        S_potential = q_over_eppmax_epp0*(-n+p-a+c+NANDNaniNcat(i));
        S_electron = g - radset*B(i)*((n*p)-(ni(i)^2)) - SRHset*(((n*p)-ni(i)^2)/((taun(i)*(p+pt(i)))+(taup(i)*(n+nt(i)))));
        S_hole     = g - radset*B(i)*((n*p)-(ni(i)^2)) - SRHset*(((n*p)-ni(i)^2)/((taun(i)*(p+pt(i)))+(taup(i)*(n+nt(i)))));
        S = [S_potential; S_electron; S_hole];
        
        if N_ionic_species == 1 || N_ionic_species == 2  % Condition for cation and anion terms
            
            F_cation = mucat(i)*(c*dVdx + kBT*(dcdx + (c*(dcdx/(DOScat(i)-c)))));

            F = [F; K_cation*mobseti*F_cation];
            
            S_cation = 0;
            S = [S; S_cation];
        end
        
        if N_ionic_species == 2     % Condition for anion terms
            
            F_anion = muani(i)*(a*-dVdx + kBT*(dadx+(a*(dadx/(DOSani(i)-a)))));
            F = [F; K_anion*mobseti*F_anion];
            
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
        elseif length(varargin) == 1 || length(varargin) >= 1 && max(max(max(varargin{1, 1}.u))) ~= 0
            switch par.N_ionic_species
                case 0
                    u0 = [interp1(icx,icsol(end,:,1),x);
                        interp1(icx,icsol(end,:,2),x);
                        interp1(icx,icsol(end,:,3),x);];
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
                if par.Rs == 0
                    Vres = 0;
                else
                    Jr = par.e*sp_r*(ur(3) - pright) - par.e*sn_r*(ur(2) - nright);
                    if par.Rs_initial
                        Vres = Jr*par.Rs*t/par.tmax;    % Initial linear sweep
                    else
                        Vres = Jr*par.Rs;
                    end
                end
                
                Pl = [-ul(1);
                    mobset*(-sn_l*(ul(2) - nleft));
                    mobset*(-sp_l*(ul(3) - pleft));];
                
                Ql = [0; 1; 1;];
                
                Pr = [-ur(1)+Vbi-Vapp+Vres;
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
