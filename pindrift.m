function solstruct = pindrift(varargin)

% A routine to test solving the diffusion and drift equations using the
% matlab pde solver. This version defines electrons as n and holes as p,
% and ions as a. V is the electric potential.
%
% Piers Barnes last modified (09/01/2016)
% Phil Calado last modified (19/11/2018)

% This version allows a previous solution to be used as the input
% conditions. If there is no input argument asssume default flat background
% condtions. If there is one argument, assume it is the previous solution
% to be used as the initial conditions. If there are two input arguments,
% assume that first are the x points from the previous solution, and the
% second is the previous solution.

if length(varargin) == 0

    p = pinParams;                         % Calls Function EParams and stores in sturcture 'params'

elseif length(varargin) == 1
    
    % Call input parameters function
    icsol = varargin{1, 1}.sol;
    icx = varargin{1, 1}.x;
    p = pinParams;                         % Calls Function EParams and stores in sturcture 'params'   

elseif length(varargin) == 2 

    if max(max(max(varargin{1, 1}.sol))) == 0

       p = varargin{2};
    
    elseif isa(varargin{2}, 'char') == 1            % Checks to see if argument is a character
        
        input_solstruct = varargin{1, 1};
        p = input_solstruct.p;
        icsol = input_solstruct.sol;
        icx = input_solstruct.x;
    
    else
        
        input_solstruct = varargin{1, 1};
        icsol = input_solstruct.sol;
        icx = input_solstruct.x;
        p = varargin{2};
    
    end

end

%% Unpack dependent properties
% Prevents recalculation of dependent properties by pdepe defined in Methods
% Can also use AbortSet in class def
        
        dcum = p.dcum;
        dEAdx = p.dEAdx;
        dIPdx = p.dIPdx;
        dN0dx = p.dN0dx;
        E0 = p.E0;
        Eg = p.Eg;
        Eif = p.Eif;
        NA = p.NA;
        ND = p.ND;
        Vbi = p.Vbi;
        n0 = p.n0;
        nleft = p.nleft;
        nright = p.nright;
        ni = p.ni; 
        p0 = p.p0;
        pleft = p.pleft;
        pright = p.pright;
        wn = p.wn;
        wp = p.wp;
        wscr = p.wscr;
        dmax = p.dcum(end);
        xx = p.xx;
        dev = p.dev;
        
%%%% Spatial mesh %%%%
if length(varargin) == 0 || max(max(max(varargin{1, 1}.sol))) == 0
    
    % Generate mesh - had problems with using the mesh generated in pc-
    % leads to very slow solving time - mesh must be being regnerated
    % constantly
    x = meshgen_x(p);
    
else
          
    x = icx;
    
end

% Define the x points to give the initial
if length(varargin) == 0 || length(varargin) == 2 && max(max(max(varargin{1, 1}.sol))) == 0
    
    icx = x;

end

p.xpoints = length(x);
%dmax = x(end);
xnm = x*1e7;        

t = meshgen_t(p);

%% Generation
% Beer Lambert
if p.OM == 1
    
    xmesh = x;  % Duplicates x the xmesh for interpolation later- there may be an easier way to do this since the beer lambert code matches the mesh.
    
    % AM15 bias light
    if p.Int ~= 0 % OM = Optical Model
        if input_solstruct.gx.AM15 == 0        % Only call BEERLAMBERT if no AM15 solution exists 
            gx.AM15 = beerlambert(p, x, 'AM15', 0, 0);
        else
            gx.AM15 = input_solstruct.gx.AM15;    % Use previous AM15 generation profile to avoid recalculation
        end
    else
        gx.AM15 = 0;
    end
    
    % laser pulse
    if p.pulseon == 1
        gx.las = beerlambert(p, x, 'laser', p.laserlambda, 0);
    else
        gx.las = 0;
    end
    
    
% Transfer Matrix
elseif p.OM == 2 && p.Int ~= 0
 
    p.genspace = x(x > p.dcum(1) & x < p.dcum(2));    % Active layer points for interpolation- this could all be implemented better but ea
    
    %% TRANSFER MATRIX NOT CURRENTLY AVAILABLE
    % Call Transfer Matrix code: [Gx1, Gx2] = TMPC1(layers, thicknesses, activeLayer1, activeLayer2)
    [Gx1S, GxLas] = TMPC1({'TiO2' 'MAPICl' 'Spiro'}, [p.d(1)+0.5*p.dint, p.d(2)+0.5*p.dint, p.d(3)+0.5*p.dint], 2, 2, p.laserlambda, p.pulsepow);
    Gx1S = Gx1S';
    GxLas = GxLas';
  
end

%% Solver options
% MaxStep = limit maximum time step size during integration
options = odeset('MaxStep', 0.1*abs(p.tmax - p.t0), 'RelTol', p.RelTol, 'AbsTol', p.AbsTol, 'NonNegative', [1,1,1,0]);

%% Call solver
% inputs with '@' are function handles to the subfunctions
% below for the: equation, initial conditions, boundary conditions
sol = pdepe(p.m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,options);

%% Subfunctions
% Set up partial differential equation (pdepe) (see MATLAB pdepe help for details of c, f, and s)
function [c,f,s,iterations] = pdex4pde(x,t,u,DuDx)

% Beer Lambert or Transfer Matrix 1 Sun
if p.OM == 1
    %AM15
    if p.Int ~= 0 
        gxAM15 = p.Int*interp1(xmesh, gx.AM15, x);
    else
        gxAM15 = 0;
    end
    %Pulse
    if p.pulseon == 1   % Applies pulse for duration set by user
        if  t >= p.pulsestart && t < p.pulselen + p.pulsestart
            gxlas = interp1(xmesh, gx.las, x);
        else
            gxlas = 0;
        end
    else
        gxlas = 0;
    end
    
    g = gxAM15 + gxlas;
        
    
elseif p.Int ~= 0 && p.OM == 2
     
      if x > dcum(1) && x <= dcum(2)
          g = p.Int*interp1(p.genspace, Gx1S, (x-dcum(1)));

      else
          g = 0;
      end
 
    % Add pulse
    if p.pulseon == 1
        if  t >= 10e-6 && t < p.pulselen + 10e-6
           if x > dcum(1) && x < dcum(2)
                lasg = p.pulseint*interp1(p.genspace, GxLas, (x-dcum(1)));
                g = g + lasg;
           end
        end
    end
  
% Uniform Generation
elseif p.OM == 0
 
           g = p.Int*p.dev.G0;

        % Add pulse
        if p.pulseon == 1
            if  t >= p.pulsestart && t < p.pulselen + p.pulsestart
                
                g = g+(p.pulseint*p.G0);
            
            end
        end
        
else
        g = 0;
        
end

% Prefactors set to 1 for time dependent components - can add other
% functions if you want to include the multiple trapping model
i = find(p.xx <= x);
i = i(end);

c = [1
     1
     1
     0];
 
if p.stats == 'Fermi'
    Dn = F.D(u(1), dev.Dnfun(i,:), dev.n_fd(i,:));
    Dp = F.D(u(2), dev.Dpfun(i,:), dev.p_fd(i,:));
elseif p.stats == 'Boltz'
    Dn = dev.mue(i)*p.kB*p.T;
    Dp = dev.muh(i)*p.kB*p.T;
end    

f = [p.mobset*(dev.mue(i)*(u(1)*(-DuDx(4)+dev.gradEA(i)-(dev.gradN0(i)*p.kB*p.T/dev.N0(i))))+(Dn*DuDx(1)));
     p.mobset*(dev.muh(i)*(u(2)*(DuDx(4)-dev.gradIP(i)-(dev.gradN0(i)*p.kB*p.T/dev.N0(i))))+(Dp*DuDx(2)));     
     p.mobseti*(dev.muion(i)*(u(3)*DuDx(4)+p.kB*p.T*(DuDx(3)+(u(3)*(DuDx(3)/(dev.DOSion(i)-u(3)))))));       % Nerst-Planck-Poisson approach ref: Borukhov 1997
     (dev.epp(i)/max(p.epp))*DuDx(4);];                                         
 
 s = [g(i) - dev.krad(i)*((u(1)*u(2))-(dev.ni(i)^2)) - p.SRHset*(((u(1)*u(2))-dev.ni(i)^2)/((dev.taun(i)*(u(2)+dev.pt(i))) + (dev.taup(i)*(u(1)+dev.nt(i)))));
      g(i) - dev.krad(i)*((u(1)*u(2))-(dev.ni(i)^2)) - p.SRHset*(((u(1)*u(2))-dev.ni(i)^2)/((dev.taun(i)*(u(2)+dev.pt(i))) + (dev.taup(i)*(u(1)+dev.nt(i)))));
      0;
      (p.q/(max(p.epp)*p.epp0))*(-u(1)+u(2)-dev.NA(i)+dev.ND(i)-dev.Nion(i)+u(3));]; 

end

% --------------------------------------------------------------------------
% Define initial conditions.
function u0 = pdex4ic(x)

if isempty(varargin) || length(varargin) >= 1 && max(max(max(varargin{1, 1}.sol))) == 0
 
i = find(p.xx <= x);
i = i(end);

u0 = [dev.n0(i);
    dev.p0(i);
    dev.Nion(i);
    dev.E0(i)];

     elseif length(varargin) == 1 || length(varargin) >= 1 && max(max(max(varargin{1, 1}.sol))) ~= 0
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

switch p.JV
    case 1
        p.Vapp = p.Vstart + ((p.Vend-p.Vstart)*t*(1/p.tmax));
    case 2 % for custom profile of voltage
        p.Vapp = p.Vapp_func(p.Vapp_params, t);
end        
    
% Open circuit condition- symmetric model
if p.OC == 1
    
    
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
    if p.BC == 0
        
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
            -ur(4) + Vbi - p.Vapp;];
        
        qr = [1;
            1;
            1;
            0];
        
    % Fixed majority charge densities at the boundaries- contact in equilibrium with etl and htl
    % Blocking electrode- zero flux for minority carriers
    elseif p.BC == 1

      
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
            -ur(4)+Vbi-p.Vapp;];
        
        qr = [0;
            1;
            1;
            0];
        
    % Non- selective contacts - fixed charge densities for majority carrier
    % and flux for minority carriers- use recombination coefficients sn_rec
    % & sp_rec to set the surface recombination velocity. pdepe
    elseif p.BC == 2
      
        pl = [-p.sn_l*(ul(1) - nleft);
            ul(2) - pleft;
            0;
            -ul(4);];
        
        ql = [1;
            0;
            1;
            0];
        
        pr = [ur(1) - nright;
            p.sp_r*(ur(2) - pright);
            0;
            -ur(4)+Vbi-p.Vapp;];
        
        qr = [0;
            1;
            1;
            0];
    
        % Flux boundary conditions for both carrier types. Leads to
        % inaccurate calculation at low currents due to the small
        % fractional change in majority carrier density at interface. May
        % be more reliable at high currents than BC2 since does not rely on
        % integrating continuity equations across the device.
    elseif p.BC == 3
      
        pl = [-p.sn_l*(ul(1) - nleft);
            -p.sp_l*(ul(2) - pleft);
            0;
            -ul(4);];
        
        ql = [1;
            1;
            1;
            0];
        
        pr = [p.sn_r*(ur(1) - nright);
            p.sp_r*(ur(2) - pright);
            0;
            -ur(4)+Vbi-p.Vapp;];
        
        qr = [1;
            1;
            1;
            0];
        
    end   
    
end

end

%% Analysis, graphing-  required to obtain J and Voc

% Readout solutions to structure
solstruct.sol = sol;            
solstruct.x = x; 
solstruct.t = t;

% include generation profile when beer lambert is included
if p.OM == 1
    solstruct.gx = gx;
end

% Store parameters structure
solstruct.p = p;

if p.OM == 2 && p.Int ~= 0
    
    solstruct.g = p.Int*interp1(p.genspace, Gx1S, (x-dcum(1)));
    
end

if p.Ana == 1
    
    [Vapp_arr, Jtotr] = pinana(solstruct, t(end)); 
       
    if p.JV == 1
        
        solstruct.Vapp = Vapp_arr;
        
    end
    
solstruct.Jtotr = Jtotr;
    
end


end
