function solstruct = df(varargin)

% Core DRIFTFUSION function- organises properties and inputs for pdepe
% A routine to test solving the diffusion and drift equations using the
% matlab pepde solver. This version defines electrons as n and holes as p,
% and ions as a. V is the electric potential.
%
% Authors: Phil Calado, Piers Barnes, Ilario Gelmetti, Ben Hillman,
% Imperial College London, 2019

% Deal with varargin arguments
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

%% Unpack dependent properties
% Prevents recalculation of dependent properties by pdepe defined in Methods
% Can also use AbortSet in class def

Vbi = par.Vbi;
nleft = par.nleft;
nright = par.nright;
pleft = par.pleft;
pright = par.pright;
xmesh = par.xx;
x_ihalf = par.x_ihalf;
devihalf = par.dev;
dev = par.dev;

%% Spatial mesh
x = xmesh;
par.xpoints = length(x);

%% Time mesh
t = meshgen_t(par);
tmesh = t;

%% Generation function
g1_fun = fun_gen(par.g1_fun_type);
g2_fun = fun_gen(par.g2_fun_type);
gxt1 = 0;
gxt2 = 0;
g = 0;

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
% Set up partial differential equation (pdepe) (see MATLAB pdepe help for details of c, f, and s)
    function [c,flux,source,iterations] = dfpde(x,t,u,DuDx)

        % Get position point
        i = find(x_ihalf <= x);
        i = i(end);

        switch par.g1_fun_type
            case 'constant'
                gxt1 = par.int1*par.gx1(i);
            otherwise
                gxt1 = g1_fun(par.g1_fun_arg, t)*par.gx1(i);
        end
        
        switch par.g2_fun_type
            case 'constant'
                gxt2 = par.int2*par.gx2(i);
            otherwise
                gxt2 = g2_fun(par.g2_fun_arg, t)*par.gx2(i);
        end

        g = gxt1 + gxt2;

        switch par.N_ionic_species
            case 1
                % Prefactors set to 1 for time dependent components - can add other
                % functions if you want to include the multiple trapping model
                c = [1;
                    1;
                    1;
                    0;];

                if par.stats == 'Fermi'
                    Dn = F.D(u(1), devihalf.Dnfun(i,:), devihalf.n_fd(i,:));
                    Dp = F.D(u(2), devihalf.Dpfun(i,:), devihalf.p_fd(i,:));
                elseif par.stats == 'Boltz'
                    Dn = devihalf.mue(i)*par.kB*par.T;
                    Dp = devihalf.muh(i)*par.kB*par.T;
                end

                flux = [par.mobset*(devihalf.mue(i)*u(1)*(-DuDx(4)+devihalf.gradEA(i))+(Dn*(DuDx(1)-((u(1)/devihalf.Nc(i))*devihalf.gradNc(i)))));
                    par.mobset*(devihalf.muh(i)*u(2)*(DuDx(4)-devihalf.gradIP(i))+(Dp*(DuDx(2)-((u(2)/devihalf.Nv(i))*devihalf.gradNv(i)))));
                    par.K_cation*par.mobseti*(devihalf.mucat(i)*(u(3)*DuDx(4)+par.kB*par.T*(DuDx(3)+(u(3)*(DuDx(3)/(devihalf.DOScat(i)-u(3)))))));       % Nerst-Planck-Poisson approach ref: Borukhov 1997
                    (devihalf.epp(i)/max(par.epp))*DuDx(4);];

                source = [g - par.radset*devihalf.krad(i)*((u(1)*u(2))-(devihalf.ni(i)^2)) - par.SRHset*(((u(1)*u(2))-devihalf.ni(i)^2)/((devihalf.taun(i)*(u(2)+devihalf.pt(i))) + (devihalf.taup(i)*(u(1)+devihalf.nt(i)))));
                    g - par.radset*devihalf.krad(i)*((u(1)*u(2))-(devihalf.ni(i)^2)) - par.SRHset*(((u(1)*u(2))-devihalf.ni(i)^2)/((devihalf.taun(i)*(u(2)+devihalf.pt(i))) + (devihalf.taup(i)*(u(1)+devihalf.nt(i)))));
                    0;
                    (par.q/(max(par.epp)*par.epp0))*(-u(1)+u(2)-devihalf.NA(i)+devihalf.ND(i)-devihalf.Nani(i)+u(3))];

            case 2
                % Prefactors set to 1 for time dependent components - can add other
                % functions if you want to include the multiple trapping model
                c = [1;
                    1;
                    1;
                    0
                    1;];

                if par.stats == 'Fermi'
                    Dn = F.D(u(1), devihalf.Dnfun(i,:), devihalf.n_fd(i,:));
                    Dp = F.D(u(2), devihalf.Dpfun(i,:), devihalf.p_fd(i,:));
                elseif par.stats == 'Boltz'
                    Dn = devihalf.mue(i)*par.kB*par.T;
                    Dp = devihalf.muh(i)*par.kB*par.T;
                end

                flux = [par.mobset*(devihalf.mue(i)*u(1)*(-DuDx(4)+devihalf.gradEA(i))+(Dn*(DuDx(1)-((u(1)/devihalf.Nc(i))*devihalf.gradNc(i)))));
                    par.mobset*(devihalf.muh(i)*u(2)*(DuDx(4)-devihalf.gradIP(i))+(Dp*(DuDx(2)-((u(2)/devihalf.Nv(i))*devihalf.gradNv(i)))));
                    par.K_cation*par.mobseti*(devihalf.mucat(i)*(u(3)*DuDx(4)+par.kB*par.T*(DuDx(3)+(u(3)*(DuDx(3)/(devihalf.DOScat(i)-u(3)))))));       % Nerst-Planck-Poisson approach ref: Borukhov 1997
                    (devihalf.epp(i)/max(par.epp))*DuDx(4);
                    par.K_anion*par.mobseti*(devihalf.muani(i)*(u(5)*-DuDx(4)+par.kB*par.T*(DuDx(5)+(u(5)*(DuDx(5)/(devihalf.DOSani(i)-u(5)))))));];

                source = [g - par.radset*devihalf.krad(i)*((u(1)*u(2))-(devihalf.ni(i)^2)) - par.SRHset*(((u(1)*u(2))-devihalf.ni(i)^2)/((devihalf.taun(i)*(u(2)+devihalf.pt(i))) + (devihalf.taup(i)*(u(1)+devihalf.nt(i)))));
                    g - par.radset*devihalf.krad(i)*((u(1)*u(2))-(devihalf.ni(i)^2)) - par.SRHset*(((u(1)*u(2))-devihalf.ni(i)^2)/((devihalf.taun(i)*(u(2)+devihalf.pt(i))) + (devihalf.taup(i)*(u(1)+devihalf.nt(i)))));
                    0;
                    (par.q/(max(par.epp)*par.epp0))*(-u(1)+u(2)-devihalf.NA(i)+devihalf.ND(i)-u(5)+u(3));
                    0;];
        end
    end

% --------------------------------------------------------------------------
% Define initial conditions.
    function u0 = dfic(x)

        if isempty(varargin) || length(varargin) >= 1 && max(max(max(varargin{1, 1}.u))) == 0

            i = find(xmesh <= x);
            i = i(end);

            switch par.N_ionic_species
                case 1
                    if length(par.dcell) == 1
                        % Single layer
                        u0 = [par.nleft*exp((x*(log(nright)-log(nleft)))/par.dcum0(end));
                            par.pleft*exp((x*(log(pright)-log(pleft)))/par.dcum0(end));
                            dev.Ncat(i);
                            (x/xmesh(end))*Vbi;];
                        % Multi-layered
                    else
                        u0 = [dev.n0(i);
                            dev.p0(i);
                            dev.Ncat(i);
                            (x/xmesh(end))*Vbi;];
                    end
                case 2
                    if length(par.dcell) == 1
                        % Single layer
                        u0 = [par.nleft*exp((x*(log(nright)-log(nleft)))/par.dcum0(end));
                            par.pleft*exp((x*(log(pright)-log(pleft)))/par.dcum0(end));
                            dev.Ncat(i);
                            (x/xmesh(end))*Vbi;
                            dev.Nani(i);];
                        % Multi-layered
                    else
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

% --------------------------------------------------------------------------

% Define boundary condtions, refer PDEPE help for the precise meaning of p
% and you l and r refer to left and right.
% in this example I am controlling the flux through the boundaries using
% the difference in concentration from equilibrium and the extraction
% coefficient.
    function [pl,ql,pr,qr] = dfbc(xl,ul,xr,ur,t)
        
        switch par.V_fun_type
            case 'constant'
                Vapp = par.Vapp;
            otherwise
                Vapp = Vapp_fun(par.V_fun_arg, t);
        end

        switch par.N_ionic_species
            case 1
                switch par.BC
                    % case 1 is obsolete
                    case 2
                        % Non- selective contacts - fixed charge densities for majority carrier
                        % and flux for minority carriers- use recombination
                        % coefficients sn_l & sp_r to set the surface recombination velocity.
                        pl = [-par.sn_l*(ul(1) - nleft);
                            ul(2) - pleft;
                            0;
                            -ul(4);];

                        ql = [1;
                            0;
                            1;
                            0;];

                        pr = [ur(1) - nright;
                            par.sp_r*(ur(2) - pright);
                            0;
                            -ur(4)+Vbi-Vapp;];

                        qr = [0;
                            1;
                            1;
                            0;];

                        % Flux boundary conditions for both carrier types. Leads to
                        % inaccurate calculation at low currents due to the small
                        % fractional change in majority carrier density at interface. May
                        % be more reliable at high currents than BC2 since does not rely on
                        % integrating continuity equations across the device.
                    case 3            
                        % Calculate series resistance voltage Vres
                        if par.Rs == 0
                            Vres = 0;
                        else
                            Jr = par.e*par.sp_r*(ur(2) - pright) - par.e*par.sn_r*(ur(1) - nright);
                            if par.Rs_initial
                                % Initial linear sweep
                                Vres = Jr*par.Rs*t/par.tmax;
                            else
                                Vres = Jr*par.Rs;
                            end
                        end

                        pl = [par.mobset*(-par.sn_l*(ul(1) - nleft));
                            par.mobset*(-par.sp_l*(ul(2) - pleft));
                            0;
                            -ul(4);];

                        ql = [1;
                            1;
                            1;
                            0;];

                        pr = [par.mobset*(par.sn_r*(ur(1) - nright));
                            par.mobset*(par.sp_r*(ur(2) - pright));
                            0;
                            -ur(4)+Vbi-Vapp+Vres;];

                        qr = [1;
                            1;
                            1;
                            0;];
                end

            case 2
                switch par.BC
                    case 2
                        % Non- selective contacts - fixed charge densities for majority carrier
                        % and flux for minority carriers- use recombination
                        % coefficients sn_l & sp_r to set the surface recombination velocity.
                        pl = [-par.sn_l*(ul(1) - nleft);
                            ul(2) - pleft;
                            0;
                            -ul(4);
                            0;];

                        ql = [1;
                            0;
                            1;
                            0;
                            1;];

                        pr = [ur(1) - nright;
                            par.sp_r*(ur(2) - pright);
                            0;
                            -ur(4)+Vbi-Vapp;
                            0;];

                        qr = [0;
                            1;
                            1;
                            0;
                            1;];

                    case 3
                        % Flux boundary conditions for both carrier types. Leads to
                        % inaccurate calculation at low currents due to the small
                        % fractional change in majority carrier density at interface. May
                        % be more reliable at high currents than BC2 since does not rely on
                        % integrating continuity equations across the device.
                        
                        % Calculate series resistance voltage Vres
                        if par.Rs == 0
                            Vres = 0;
                        else
                            Jr = par.e*par.sp_r*(ur(2) - pright) - par.e*par.sn_r*(ur(1) - nright);
                            if par.Rs_initial
                                % Initial linear sweep
                                Vres = Jr*par.Rs*t/par.tmax;
                            else
                                Vres = Jr*par.Rs;
                            end
                        end

                        pl = [par.mobset*(-par.sn_l*(ul(1) - nleft));
                            par.mobset*(-par.sp_l*(ul(2) - pleft));
                            0;
                            -ul(4);
                            0;];

                        ql = [1;
                            1;
                            1;
                            0;
                            1;];

                        pr = [par.mobset*(par.sn_r*(ur(1) - nright));
                            par.mobset*(par.sp_r*(ur(2) - pright));
                            0;
                            -ur(4)+Vbi-Vapp+Vres;
                            0;];

                        qr = [1;
                            1;
                            1;
                            0;
                            1;];
                end
        end
    end

% Store final voltage reading
par.Vapp = Vapp;

% Readout solutions to structure
solstruct.u = u;
solstruct.x = x;
solstruct.t = t;

% Store parameters structure
solstruct.par = par;

end
