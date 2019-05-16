function solstruct = df(varargin)

% Core DRIFTFUSION function- organises properties and inputs for pdepe
% A routine to test solving the diffusion and drift equations using the
% matlab pepde solver. This version defines electrons as n and holes as p,
% and ions as a. V is the electric potential.
%
% Authors: Phil Calado, Piers Barnes, Ilario Gelmetti, Ben Hillman,
% Imperial College London, 2018

% Variable input arguments

% With no input argument, parameters
if length(varargin) == 0
    
    par = pc;
    
elseif length(varargin) == 1
    
    % Call input parameters function
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

dcum = par.dcum;
E0 = par.E0;
Eg = par.Eg;
Eif = par.Eif;
NA = par.NA;
ND = par.ND;
Vbi = par.Vbi;
n0 = par.n0;
nleft = par.nleft;
nright = par.nright;
ni = par.ni;
p0 = par.p0;
pleft = par.pleft;
pright = par.pright;
dmax = par.dcum(end);
xx = par.xx;
dev = par.dev;

%%%% Spatial mesh %%%%
if length(varargin) == 0 || max(max(max(varargin{1, 1}.u))) == 0
    
    % Generate mesh - had problems with using the mesh generated in pc-
    % leads to very slow solving time - mesh must be being regnerated
    % constantly
    x = meshgen_x(par);
    
else
    
    x = icx;
    
end

% Define the x points to give the initial
if length(varargin) == 0 || length(varargin) == 2 && max(max(max(varargin{1, 1}.u))) == 0
    
    icx = x;
    
end

par.xpoints = length(x);
%dmax = x(end);
xnm = x*1e7;

t = meshgen_t(par);

%% Generation
% Beer Lambert
if par.OM == 1
        
    % AM15 bias light
    if par.Int ~= 0 % OM = Optical Model
        if input_solstruct.gx.AM15 == 0        % Only call BEERLAMBERT if no AM15 solution exists
            gx.AM15 = beerlambert(par, x, 'AM15', 0, 0);
        else
            gx.AM15 = input_solstruct.gx.AM15;    % Use previous AM15 generation profile to avoid recalculation
        end
    else
        gx.AM15 = zeros(length(x), 1);
    end
    
    % laser pulse
    if par.pulseon == 1
        gx.las = beerlambert(par, x, 'laser', par.laserlambda, 0);
    else
        gx.las = zeros(length(x), 1);
    end
    
end

% Transfer Matrix
if par.OM == 2 && par.Int ~= 0
    
    par.genspace = x(x > par.dcum(1) & x < par.dcum(2));    % Active layer points for interpolation- this could all be implemented better but ea
    
    %% TRANSFER MATRIX NOT CURRENTLY AVAILABLE
    % Call Transfer Matrix code: [Gx1, Gx2] = TMPC1(layers, thicknesses, activeLayer1, activeLayer2)
    [Gx1S, GxLas] = TMPC1({'TiO2' 'MAPICl' 'Spiro'}, [par.d(1)+0.5*par.dint, par.d(2)+0.5*par.dint, par.d(3)+0.5*par.dint], 2, 2, par.laserlambda, par.pulsepow);
    Gx1S = Gx1S';
    GxLas = GxLas';
    
end

%% Solver options
% MaxStep = limit maximum time step size during integration
options = odeset('MaxStep', 0.1*abs(par.tmax - par.t0), 'RelTol', par.RelTol, 'AbsTol', par.AbsTol, 'NonNegative', [1,1,1,0]);

%% Call solver
% inputs with '@' are function handles to the subfunctions
% below for the: equation, initial conditions, boundary conditions
u = pdepe(par.m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,options);

%% Subfunctions
% Set up partial differential equation (pdepe) (see MATLAB pdepe help for details of c, f, and s)
    function [c,f,s,iterations] = pdex4pde(x,t,u,DuDx)
        
        i = find(par.xx <= x);
        i = i(end);
        
        % Uniform Generation
        switch par.OM
            
            case 0
                if par.Int ~= 0
                    g = par.Int*dev.G0(i);
                else
                    g = 0;
                end
                % Add pulse
                if par.pulseon == 1
                    if  t >= par.pulsestart && t < par.pulselen + par.pulsestart
                        
                        g = g+(par.pulseint.*dev.G0);
                        
                    end
                end
                
            % Beer Lambert
            case 1
                %AM15
                    gxAM15 = par.Int*gx.AM15(i);
                %Pulse
                if par.pulseon == 1   % Applies pulse for duration set by user
                    if  t >= par.pulsestart && t < par.pulselen + par.pulsestart
                        gxlas = gx.las(i);
                    else
                        gxlas = 0;
                    end
                else
                    gxlas = 0;
                end
                
                g = gxAM15 + gxlas;
                
            % Transfer Matrix- not currently working!
            case 2
                
                if par.Int ~= 0
                    
                    if x > dcum(1) && x <= dcum(2)
                        g = par.Int*interp1(par.genspace, Gx1S, (x-dcum(1)));
                    end
                else
                    g = 0;
                end
                
                % Add pulse
                if par.pulseon == 1
                    if  t >= 10e-6 && t < par.pulselen + 10e-6
                        if x > dcum(1) && x < dcum(2)
                            lasg = par.pulseint*interp1(par.genspace, GxLas, (x-dcum(1)));
                            g = g + lasg;
                        end
                    end
                end
                
        end
        
        % Prefactors set to 1 for time dependent components - can add other
        % functions if you want to include the multiple trapping model
        c = [1;
            1;
            1;
            0];
        
        if par.stats == 'Fermi'
            Dn = F.D(u(1), dev.Dnfun(i,:), dev.n_fd(i,:));
            Dp = F.D(u(2), dev.Dpfun(i,:), dev.p_fd(i,:));
        elseif par.stats == 'Boltz'
            Dn = dev.mue(i)*par.kB*par.T;
            Dp = dev.muh(i)*par.kB*par.T;
        end
        
        f = [par.mobset*(dev.mue(i)*u(1)*(-DuDx(4)+dev.gradEA(i))+(Dn*(DuDx(1)-((u(1)/dev.Nc(i))*dev.gradNc(i)))));
            par.mobset*(dev.muh(i)*u(2)*(DuDx(4)-dev.gradIP(i))+(Dp*(DuDx(2)-((u(2)/dev.Nv(i))*dev.gradNv(i)))));
            par.mobseti*(dev.muion(i)*(u(3)*DuDx(4)+par.kB*par.T*(DuDx(3))));%+(u(3)*(DuDx(3)/(dev.DOSion(i)-u(3)))))));       % Nerst-Planck-Poisson approach ref: Borukhov 1997
            (dev.epp(i)/max(par.epp))*DuDx(4);];
        
        s = [g - dev.krad(i)*((u(1)*u(2))-(dev.ni(i)^2)) - par.SRHset*(((u(1)*u(2))-dev.ni(i)^2)/((dev.taun(i)*(u(2)+dev.pt(i))) + (dev.taup(i)*(u(1)+dev.nt(i)))));
            g - dev.krad(i)*((u(1)*u(2))-(dev.ni(i)^2)) - par.SRHset*(((u(1)*u(2))-dev.ni(i)^2)/((dev.taun(i)*(u(2)+dev.pt(i))) + (dev.taup(i)*(u(1)+dev.nt(i)))));
            0;
            (par.q/(max(par.epp)*par.epp0))*(-u(1)+u(2)-dev.NA(i)+dev.ND(i)-dev.Nion(i)+u(3));];
        
    end

% --------------------------------------------------------------------------
% Define initial conditions.
    function u0 = pdex4ic(x)
        
        if isempty(varargin) || length(varargin) >= 1 && max(max(max(varargin{1, 1}.u))) == 0
            
            i = find(par.xx <= x);
            i = i(end);
            
            u0 = [dev.n0(i);
                dev.p0(i);
                dev.Nion(i);
                dev.E0(i)];
            
        elseif length(varargin) == 1 || length(varargin) >= 1 && max(max(max(varargin{1, 1}.u))) ~= 0
            % insert previous solution and interpolate the x points
            u0 = [interp1(icx,icsol(end,:,1),x)
                interp1(icx,icsol(end,:,2),x)
                interp1(icx,icsol(end,:,3),x)
                interp1(icx,icsol(end,:,4),x)];
            
        end
        
    end

% --------------------------------------------------------------------------

% Define boundary condtions, refer pdepe help for the precise meaning of p
% and you l and r refer to left and right.
% in this example I am controlling the flux through the boundaries using
% the difference in concentration from equilibrium and the extraction
% coefficient.
    function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
        
        switch par.JV
            case 1
                par.Vapp = par.Vstart + ((par.Vend-par.Vstart)*t*(1/par.tmax));
            case 2 % for custom profile of voltage
                par.Vapp = par.Vapp_func(par.Vapp_params, t);
        end
        
        % Open circuit condition- symmetric model
        if par.OC == 1
            
            
            pl = [0;
                0;
                0;
                -ul(4)];
            
            ql = [1;
                1;
                1;
                0];
            
            pr = [0;
                0;
                0;
                -ur(4)];
            
            qr = [1;
                1;
                1;
                0];
            
        else
            %% Closed circuit condition
            
            % Zero current - rarely used but can be useful to switch off currents before switching
            % to OC in procedures.
            if par.BC == 0
                
                pl = [0;
                    0;
                    0;
                    -ul(4)];
                
                ql = [1;
                    1;
                    1;
                    0];
                
                pr = [0;
                    0;
                    0;
                    -ur(4) + Vbi - par.Vapp;];
                
                qr = [1;
                    1;
                    1;
                    0];
                
                % Fixed majority charge densities at the boundaries- contact in equilibrium with etl and htl
                % Blocking electrode- zero flux for minority carriers
            elseif par.BC == 1
                
                
                pl = [0;
                    (ul(2)-pleft);
                    0;
                    -ul(4);];
                
                ql = [1;
                    0;
                    1;
                    0];
                
                pr = [(ur(1)-nright);
                    0;
                    0;
                    -ur(4)+Vbi-par.Vapp;];
                
                qr = [0;
                    1;
                    1;
                    0];
                
                % Non- selective contacts - fixed charge densities for majority carrier
                % and flux for minority carriers- use recombination coefficients sn_rec
                % & sp_rec to set the surface recombination velocity. pdepe
            elseif par.BC == 2
                
                pl = [-par.sn_l*(ul(1) - nleft);
                    ul(2) - pleft;
                    0;
                    -ul(4);];
                
                ql = [1;
                    0;
                    1;
                    0];
                
                pr = [ur(1) - nright;
                    par.sp_r*(ur(2) - pright);
                    0;
                    -ur(4)+Vbi-par.Vapp;];
                
                qr = [0;
                    1;
                    1;
                    0];
                
                % Flux boundary conditions for both carrier types. Leads to
                % inaccurate calculation at low currents due to the small
                % fractional change in majority carrier density at interface. May
                % be more reliable at high currents than BC2 since does not rely on
                % integrating continuity equations across the device.
            elseif par.BC == 3
                Jn_r = -par.e*par.sn_r*(ur(1) - nright);
                Jp_r = par.e*par.sp_r*(ur(2) - pright);
                
                Jr = Jn_r+Jp_r;
                Vres = Jr*par.Rs;
                
                pl = [-par.sn_l*(ul(1) - nleft);
                    -par.sp_l*(ul(2) - pleft);
                    0;
                    -ul(4);];
                
                ql = [1;
                    1;
                    1;
                    0];
                
                pr = [par.sn_r*(ur(1) - nright);
                    par.sp_r*(ur(2) - pright);
                    0;
                    -ur(4)+Vbi-par.Vapp+Vres;];
                
                qr = [1;
                    1;
                    1;
                    0];
                
            end
            
        end
        
    end

%% Analysis, graphing-  required to obtain J and Voc

% Readout solutions to structure
solstruct.u = u;
solstruct.x = x;
solstruct.t = t;

% include generation profile when beer lambert is included
if par.OM == 1
    solstruct.gx = gx;
end

% Store parameters structure
solstruct.par = par;

if par.OM == 2 && par.Int ~= 0
    solstruct.g = par.Int*interp1(par.genspace, Gx1S, (x-dcum(1)));
end

end
