function solstruct = pindrift(varargin);

% Requires v2struct toolbox for unpacking parameters structure
% IMPORTANT! Currently uses parameters from pinParams - all variables must be
% declared in pinDrift (line 52)

% A routine to test solving the diffusion and drift equations using the
% matlab pde solver. This version defines electrons as n and holes as p,
% and ions as a. V is the electric potential.
%
% Piers Barnes last modified (09/01/2016)
% Phil Calado last modified (07/07/2016)

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

        Bn = p.Bn;
        Bp = p.Bp;
        E0 = p.E0;
        Eg = p.Eg;
        Eif = p.Eif;
        Et = p.Et;
        NA = p.NA;
        ND = p.ND;
        Vbi = p.Vbi;
        n0 = p.n0;
        ni = p.ni;
        nt = p.nt;     
        p0 = p.p0;
        pt = p.pt;
        wn = p.wn;
        wp = p.wp;
        wscr = p.wscr;

%%%% Spatial mesh %%%%
if length(varargin) == 0 || max(max(max(varargin{1, 1}.sol))) == 0
    
    % Edit meshes in mesh gen
    x = meshgen_x(p);
    
        if p.OC == 1
        
        % Mirror the mesh for symmetric model - symmetry point is an additional
        % point at device length + 1e-7
        x1 = x;
        x2 = p.xmax - fliplr(x) + x(end);
        x2 = x2(2:end);                 % Delete initial point to ensure symmetry
        x = [x1, x2];
        
        end

else
          
    x = icx;
    
end

% Define the x points to give the initial
if length(varargin) == 0 || length(varargin) == 2 && max(max(max(varargin{1, 1}.sol))) == 0
    
    icx = x;

end

p.xpoints = length(x);
p.xmax = x(end);
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
    [Gx1S, GxLas] = TMPC1({'TiO2' 'MAPICl' 'Spiro'}, [p.d(1) p.d(2) p.d(3)], 2, 2, p.laserlambda, p.pulsepow);
    Gx1S = Gx1S';
    GxLas = GxLas';
  
end

%% Solver options
% MaxStep = limit maximum time step size during integration
options = odeset('MaxStep', 0.1*abs(p.tmax - p.t0), 'RelTol', p.RelTol, 'AbsTol', p.AbsTol);

%% Call solver
% inputs with '@' are function handles to the subfunctions
% below for the: equation, initial conditions, boundary conditions
sol = pdepe(p.m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,options);

%% Subfunctions
% Set up partial differential equation (pdepe) (see MATLAB pdepe help for details of c, f, and s)
function [c,f,s,iterations] = pdex4pde(x,t,u,DuDx)

% Open circuit condition- symmetric model
if (p.OC ==1)
    
    if x > p.xmax/2

        x = p.xmax - x;
        gradcoeff = -1;
       
    else
        
        gradcoeff = 1;
        
    end

else
    
    gradcoeff = 1;
    
end

% Beer Lambert or Transfer Matrix 1 Sun
%% NOT CURRENTLY WORKING
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
     
      if x > p.tp && x < (p.tp+p.ti) 
          g = p.Int*interp1(p.genspace, Gx1S, (x-p.tp));

      else
          g = 0;
      end
 
    % Add pulse
    if p.pulseon == 1
        if  t >= 10e-6 && t < p.pulselen + 10e-6
           if x > p.tp && x < (p.tp+p.ti)
                lasg = p.pulseint*interp1(p.genspace, GxLas, (x-p.tp));
                g = g + lasg;
           end
        end
    end
  
% Uniform Generation
elseif p.OM == 0
      
      if p.Int ~= 0 && x > p.tp && x < (p.tp+p.ti)    
           g = p.Int*p.G0;
      else
           g = 0;
      end
        
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
c = [1
     1
     1
     0];

% p-type
if x <= p.dcum(1)
 
   f = [(Bn(1)*p.mue(1)*(u(1)*-DuDx(4)+p.kB*p.T*DuDx(1)));
     (Bp(1)*p.muh(1)*(u(2)*DuDx(4)+p.kB*p.T*DuDx(2)));     
     0;
     p.epp(1)*DuDx(4);];                                  

 s = [ - p.kradp*((u(1)*Bn(1)*u(2)*Bp(1))-(ni(1)^2));% - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+ptp)) + (taup(1)*(u(1)+ntp)))); %; %- klincon*min((u(1)- htln0), (u(2)- htlp0)); % 
       - p.kradp*((u(1)*Bn(1)*u(2)*Bp(1))-(ni(1)^2));% - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+ptp)) + (taup(1)*(u(1)+ntp)))); %- kradhtl*((u(1)*u(2))-(ni^2)); %- klincon*min((u(1)- htln0), (u(2)- htlp0)); % - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+pthtl)) + (taup(1)*(u(1)+nthtl))));
      0;
      (p.q/p.epp0)*(-(u(1)*Bn(1))+(u(2)*Bp(1))-NA(1)+ND(1)-p.NI+u(3));];

% Intrinsic Recombination zone
elseif x >  p.dcum(1) && x <= p.dcum(1) +p.tint
% Virtual charge carrier densities based on Fermi level equilibrium
f = [(Bn(2)*p.mue(2)*((u(1)*-DuDx(4))+(p.kB*p.T*DuDx(1))));
     (Bp(2)*p.muh(2)*((u(2)*DuDx(4))+(p.kB*p.T*DuDx(2))));     
     (p.mui*(u(3)*DuDx(4)+p.kB*p.T*DuDx(3)));
     p.epp(2)*DuDx(4);];                                      

 s = [g - p.kradi*((u(1)*Bn(2)*u(2)*Bp(2))-(ni(2)^2)) - (((u(1)*Bn(2)*u(2)*Bp(1))-ni(2)^2)/((p.taun(1)*(u(2)*Bp(1)+pt(2))) + (p.taup(1)*(u(1)*Bn(2)+nt(2))))); %- krad*((u(1)*u(2))-(ni^2));  % - klin*min((u(1)- ni), (u(2)- ni)); % - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+ptrap)) + (taup(1)*(u(1)+ntrap))));
      g - p.kradi*((u(1)*Bn(2)*u(2)*Bp(2))-(ni(2)^2)) - (((u(1)*Bn(2)*u(2)*Bp(1))-ni(2)^2)/((p.taun(1)*(u(2)*Bp(1)+pt(2))) + (p.taup(1)*(u(1)*Bn(2)+nt(2))))); %- krad*((u(1)*u(2))-(ni^2));  % - klin*min((u(1)- ni), (u(2)- ni)); % - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+ptrap)) + (taup(1)*(u(1)+ntrap))));
      0;
      (p.q/p.epp0)*(-(u(1)*Bn(2))+(u(2)*Bp(2))-p.NI+u(3));]; 

% Intrinsic
elseif x >  p.dcum(1) +p.tint && x < p.dcum(2) - p.tint
% Virtual charge carrier densities based on Fermi level equilibrium
f = [(Bn(2)*p.mue(2)*((u(1)*-DuDx(4))+(p.kB*p.T*DuDx(1))));
     (Bp(2)*p.muh(2)*((u(2)*DuDx(4))+(p.kB*p.T*DuDx(2))));     
     (p.mui*(u(3)*DuDx(4)+p.kB*p.T*DuDx(3)));
     p.epp(2)*DuDx(4);];                                      

 s = [g - p.kradi*((u(1)*Bn(2)*u(2)*Bp(2))-(ni(2)^2));% - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+pthtl)) + (taup(1)*(u(1)+nthtl)))); %- krad*((u(1)*u(2))-(ni^2));  % - klin*min((u(1)- ni), (u(2)- ni)); % - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+ptrap)) + (taup(1)*(u(1)+ntrap))));
      g - p.kradi*((u(1)*Bn(2)*u(2)*Bp(2))-(ni(2)^2));% - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+pthtl)) + (taup(1)*(u(1)+nthtl)))); %- krad*((u(1)*u(2))-(ni^2));  % - klin*min((u(1)- ni), (u(2)- ni)); % - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+ptrap)) + (taup(1)*(u(1)+ntrap))));
      0;
      (p.q/p.epp0)*(-(u(1)*Bn(2))+(u(2)*Bp(2))-p.NI+u(3));]; 
 
% Intrinsic Recombination zone
elseif x >=  p.dcum(2) - p.tint && x < p.dcum(2) 
% Virtual charge carrier densities based on Fermi level equilibrium
f = [(Bn(2)*p.mue(2)*((u(1)*-DuDx(4))+(p.kB*p.T*DuDx(1))));
     (Bp(2)*p.muh(2)*((u(2)*DuDx(4))+(p.kB*p.T*DuDx(2))));     
     (p.mui*(u(3)*DuDx(4)+p.kB*p.T*DuDx(3)));
     p.epp(2)*DuDx(4);];                                      

 s = [g - p.kradi*((u(1)*Bn(2)*u(2)*Bp(2))-(ni(2)^2)) - (((u(1)*Bn(3)*u(2)*Bp(2))-ni(2)^2)/((p.taun(3)*(u(2)*Bp(2)+pt(2))) + (p.taup(3)*(u(1)*Bn(3)+nt(2))))); %- krad*((u(1)*u(2))-(ni^2));  % - klin*min((u(1)- ni), (u(2)- ni)); % - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+ptrap)) + (taup(1)*(u(1)+ntrap))));
      g - p.kradi*((u(1)*Bn(2)*u(2)*Bp(2))-(ni(2)^2)) - (((u(1)*Bn(3)*u(2)*Bp(2))-ni(2)^2)/((p.taun(3)*(u(2)*Bp(2)+pt(2))) + (p.taup(3)*(u(1)*Bn(3)+nt(2))))); %- krad*((u(1)*u(2))-(ni^2));  % - klin*min((u(1)- ni), (u(2)- ni)); % - (((u(1)*u(2))-ni^2)/((taun(1)*(u(2)+ptrap)) + (taup(1)*(u(1)+ntrap))));
      0;
      (p.q/p.epp0)*(-(u(1)*Bn(2))+(u(2)*Bp(2))-p.NI+u(3));]; 
  
% n-type
elseif x >= p.dcum(2) &&  x <= p.xmax
    
f = [(Bn(3)*p.mue(3)*((u(1)*-DuDx(4))+(p.kB*p.T*DuDx(1))));
     (Bp(3)*p.muh(3)*((u(2)*DuDx(4))+(p.kB*p.T*DuDx(2))));      
     0;
     p.epp(3)*DuDx(4)];                                      

s = [ - p.kradn*((u(1)*Bn(3)*u(2)*Bp(3))-(ni(3)^2));% - (((u(1)*Bn(3)*u(2)*Bp(3))-ni(3)^2)/((p.taun(3)*(u(2)*Bp(3)+ptn)) + (taup(3)*(u(1)*Bn(3)+ntn))));  %- kradetl*((u(1)*u(2))-(ni^2)); %- klincon*min((u(1)- etln0), (u(2)- etlp0)); %  - (((u(1)*u(2))-ni^2)/((taun(3)*(u(2)+ptetl)) + (taup(3)*(u(1)+ntetl))));
      - p.kradn*((u(1)*Bn(3)*u(2)*Bp(3))-(ni(3)^2));% - (((u(1)*Bn(3)*u(2)*Bp(3))-ni(3)^2)/((p.taun(1)*(u(2)*Bp(3)+ptn)) + (taup(1)*(u(1)*Bn(3)+ntn))));   %- kradetl*((u(1)*u(2))-(ni^2)); % - klincon*min((u(1)- etln0), (u(2)- etlp0)); %- (((u(1)*u(2))-ni^2)/((taun(3)*(u(2)+ptetl)) + (taup(3)*(u(1)+ntetl))));
      0;
      (p.q/p.epp0)*(-(u(1)*Bn(3))+(u(2)*Bp(3))-p.NI+u(3)+ND(3)-NA(3));];

end

end

% --------------------------------------------------------------------------

% Define initial conditions.
function u0 = pdex4ic(x)

% Open circuit condition- symmetric model
if (p.OC ==1)
    
    if x > p.xmax/2

        x = p.xmax - x;

    end
    
end

if length(varargin) == 0 || length(varargin) >= 1 && max(max(max(varargin{1, 1}.sol))) == 0
 
    
%         % Intrinsic
%     u0  = [ni;
%          ni;
%          NI;
%          Ei]; 
%      
%     % p-type
%     if x < (tp+tscr - wp)
%     
%        u0 =  [n0(1);
%              p0(1);
%              0;
%              0];  
% 
%     % p-type SCR    
%     elseif  x >= (tp+tscr - wp) && x < (tp+tscr)
% 
%         u0 = [N0(1)*exp(q*(Efnside + EAp + q*((((q*NA)/(2*p.epp(2)*p.epp0))*(x-tp-tscr+wp)^2)))/(kB*T));                            %ni*exp((Efnside - (-q*((((q*NA)/(2*epp(1)))*(x-tp+wp)^2))))/(kB*T));
%               N0(1)*exp(-q*(q*((((q*NA)/(2*p.epp(2)*p.epp0))*(x-tp-tscr+wp)^2)) + EAp + Egp + Efpside)/(kB*T)) ;
%               0;
%               (((q*NA)/(2*p.epp(2)*p.epp0))*(x-tp-tscr+wp)^2)];
% 
%     % Intrinsic
% 
%     elseif x >= tp+tscr && x < tp+ ti+tscr
%         u0 =  [N0(2)*exp(q*(Efnside + EAp + q*(((x - tp-tscr)*((1/ti)*(Vbi - ((q*NA*wp^2)/(2*p.epp(2)*p.epp0)) - ((q*ND*wn^2)/(2*p.epp(2)*p.epp0))))) + ((q*NA*wp^2)/(2*p.epp(2)*p.epp0))))/(kB*T));
%                 N0(2)*exp(-q*(q*(((x - tp-tscr)*((1/ti)*(Vbi - ((q*NA*wp^2)/(2*p.epp(2)*p.epp0)) - ((q*ND*wn^2)/(2*p.epp(2)*p.epp0))))) + ((q*NA*wp^2)/(2*p.epp(2)*p.epp0))) + EAp + Egp + Efpside)/(kB*T)) ;
%                NI;
%                 ((x - tp-tscr)*((1/ti)*(Vbi - ((q*NA*wp^2)/(2*p.epp(2)*p.epp0)) - ((q*ND*wn^2)/(2*p.epp(2)*p.epp0))))) + ((q*NA*wp^2)/(2*p.epp(2)*p.epp0)) ;];
% 
% %      % n-type SCR    
%     elseif  x >= (tp+tscr+ti) && x < (tp +tscr + wn +ti)
% 
%         u0 = [N0(3)*exp(q*(Efnside + EAn*((x-tscr-ti-tp)^2)/wn^2 + EAp*((x-tscr-ti-tp-wn)^2)+ q*(((-(q*ND)/(2*p.epp(2)*p.epp0)*(x-tscr-ti-tp-wn)^2 )+ Vbi)))/(kB*T));
%               N0(3)*exp(-q*(q*((((-(q*ND)/(2*p.epp(2)*p.epp0))*(x-tscr-ti-tp-wn)^2) + Vbi)) + EAn*(x-tscr-ti-tp)^2/wn^2 + EAp*(x-tscr-ti-tp-wn)^2+ + Egp + Efpside)/(kB*T));
%               0; 
%               (((-(q*ND)/(2*p.epp(2)*p.epp0))*(x-tp-ti - tscr -wn)^2) + Vbi)];
% 
%     elseif x >= (tp + tscr + wn +ti )
% 
%          u0 = [n0(3);
%                p0(3);
%                0;
%                Vbi];
% 
%     end

% 
%      u0  = [ni(2)/(p.N0(2)/p.N0(1))*exp((p.EA(1)-p.EA(2))/(p.kB*p.T));
%             ni(2)/(p.N0(2)/p.N0(1))*exp((p.IP(2)-p.IP(1))/(p.kB*p.T));
%             p.NI;
%             Eif(2)];
        
%     elseif x >=  p.tp + p.ti 
%      
%               
%     u0  = [n0(3)/(p.N0(3)/p.N0(1))*exp((p.EA(1)-p.EA(3))/(p.kB*p.T));
%            p0(3)/(p.N0(3)/p.N0(1))*exp((p.IP(3)-p.IP(1))/(p.kB*p.T));

    if x <=  p.tp
    
     u0  = [n0(1)/Bn(1);
            p0(1)/Bp(1);
            p.NI;
            p.E0(1)];
        
    elseif x > p.tp  && x < p.tp + p.ti

     u0  = [ni(2)/Bn(2);
            ni(2)/Bp(2);
            p.NI;
            Eif(2)];
        
    elseif x >=  p.tp + p.ti
                   
    u0  = [n0(3)/Bn(3);
           p0(3)/Bp(3);
           p.NI;
           p.E0(3)];
    
    end


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

if p.JV == 1
        
    p.Vapp = p.Vstart + ((p.Vend-p.Vstart)*t*(1/p.tmax));
    
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
            (ul(2)*Bp(1)-p0(1));
            0;
            -ul(4);];
        
        ql = [1;
            0;
            1;
            0];
        
        pr = [(ur(1)*Bn(3)-n0(3));
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
      
        pl = [-p.sn_l*(ul(1)*Bn(1) - n0(1));
            ul(2)*Bp(1) - p0(1);
            0;
            -ul(4);];
        
        ql = [1;
            0;
            1;
            0];
        
        pr = [ur(1)*Bn(3) - n0(3);
            p.sp_r*(ur(2)*Bp(3) - p0(3));
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
      
        pl = [-p.sn_l*(ul(1)*Bn(1) - n0(1));
            -p.sp_l*(ul(2)*Bp(1) - p0(1));
            0;
            -ul(4);];
        
        ql = [1;
            1;
            1;
            0];
        
        pr = [p.sn_r*(ur(1)*Bn(3) - n0(3));
            p.sp_r*(ur(2)*Bp(3) - p0(3));
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
    
    solstruct.g = p.Int*interp1(p.genspace, Gx1S, (x-p.tp));
    
end

if p.Ana == 1
    
    [Voc, Vapp_arr, Jtotr] = pinAna(solstruct, t(end));
    
    if p.OC == 1
        
        solstruct.Voc = Voc;
        
    end
    
    if p.JV == 1
        
        solstruct.Vapp = Vapp_arr;
        
    end
    
    solstruct.Jtotr = Jtotr;
    
end


end
