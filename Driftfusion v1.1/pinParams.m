function [params] = pinParams
% Generates parameters structure for pinDrift

% Physical constants
p.kB = 8.6173324e-5;    % Boltzmann constant [eV K^-1]
p.T = 300;              % Temperature [K]
p.epp0 = 552434;        % [e^2 eV^-1 cm^-1] -Checked (02-11-15)
p.q = 1;                % in e
p.e = 1.61917e-19;      % Charge of an electron in Coulombs for current calculations

% Device Dimensions [cm]
p.tp = 200e-7;         % p-type layer thickness
p.pp = 60;             % p-type layer points
p.ti = 400e-7;         % Intrinsic layer thickness
p.pii = 100;           % Intrinsic points
p.tn = 200e-7;         % n-type thickness
p.pn = 60;             % n-type points
p.tint = 20e-7;      % 0.5x Interfacial region thickness (x_mesh_type = 3)
p.pint = 40;         % 0.5x Interfacial points (x_mesh_type = 3)
p.tscr = 50e-7;
p.pscr = 50;
p.pepe = 20;           % electrode interface points
p.te = 10e-7;          % electrode interface thickness
p.deltax = p.tint/p.pint;    % spacing in the interfacial region- requires for mesh generation

% Parameters for spatial mesh of solution points - see meshgen_x for
% xmesh_type specification
p.xmesh_type = 4; 
p.xmax = p.tp + p.ti + p.tn;      % cm

if p.xmesh_type == 1 || p.xmesh_type == 5

    p.x0 = 0;

else
    
    p.x0 = p.xmax/1e3;

end

% General Parameters
p.OC = 0;                 % Closed circuit = 0, Open Circuit = 1 
p.Int = 1;                % Bias Light intensity (Suns Eq.)
p.G0 = 2.5e21;            % Uniform generation rate @ 1 Sun
p.tmax = 1e-3;            % Time
p.pulseon = 0;            % Switch pulse on TPC or TPV
p.Vapp = 0;               % Applied bias
p.BC = 1;                 % Boundary Conditions. Must be set to one for first solution
p.figson = 1;             % Toggle figures on/off
p.meshx_figon = 0;        % Toggles x-mesh figures on/off
p.mesht_figon = 0;        % Toggles t-mesh figures on/off
p.side = 1;               % illumination side 1 = EE, 2 = SE
p.calcJ = 0;              % Calculates Currents- slows down solving calcJ = 1, calculates DD currents at every position, calcJ = 2, calculates DD at boundary.
p.mobset = 1;             % Switch on/off electron hole mobility- MUST BE SET TO ZERO FOR INITIAL SOLUTION
p.JV = 0;                 % Toggle run JV scan on/off
p.Ana = 1;                % Toggle on/off analysis

% OM = Optical Model
% Current only uniform generation functionality is avaiable- this will be
% updated in future versions.
% 0 = Uniform Generation 
p.OM  = 0;

% Parameters for time mesh
p.tpoints = 200;              % Number of time points
p.tmesh_type = 2;             % Mesh type- for use with meshgen_t
p.t0 = p.tmax/1e4;
p.deltat = p.tmax/(1e2*p.tpoints);

%%%%%%%%%% CURRENT VOLTAGE SCAN %%%%%%%%%%%%
% NOT 100% reliable- requires further development and investigation
p.Vstart = 0;
p.Vend = 1.2;
p.JVscan_rate = 1;        % JV scan rate (Vs-1)
p.JVscan_pnts = 600;

%%%%%%%% TRANSIENT SETTINGS %%%%%%%%
if p.pulseon && p.OC       % Record length for TPV
    
    p.pulselen = 1e-6;       % Transient pulse length
    p.pulseint = 10;          % Transient pulse intensity- for BL and TM models, 100 mW/cm2 assumed
    p.pulsestart = 1e-7;     % Time before pulse- required for baseline
    p.tmesh_type = 3;
    p.tmax = 100e-6;
    p.t0 = tmax/1e6;
    p.tpoints = 400;
    p.deltat = tmax/(1e6*tpoints);
    p.JV = 0;
    p.calcJ = 0;

end    

if p.pulseon && ~ p.OC   % Record length for TPC
    
    p.pulselen = 1e-6;       % Transient pulse length
    p.pulseint = 0.2;            % Transient pulse intensity- for BL and TM models, 100 mW/cm2 assumed
    p.pulsestart = 1e-7;
    p.tmesh_type = 3;
    p.tmax = 1e-6;  
    p.t0 = 0;%tmax/1000;
    p.tpoints = 400;
    p.deltat = p.tmax/(1e4*p.tpoints);
    p.calcJ = 2;
    p.JV = 0;
    
end

%%%%%%%%%%% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%
if ~ p.mobset
   
    p.mue_i = 0;     % [Vcm-2s-1] electron mobility
    p.muh_i = 0;     % hole mobility
    p.mue_p = 0;
    p.muh_p = 0;
    p.mue_n = 0;
    p.muh_n = 0;
    
else
    
    p.mue_i = 20;          % electron mobility
    p.muh_i = 20;      % hole mobility
    p.mue_p = p.mue_i;
    p.muh_p = p.mue_i;
    p.mue_n = p.mue_i;
    p.muh_n = p.mue_i;
    
end

p.mui = 1e-10; % ion mobility

p.eppp = 20*p.epp0;         % Dielectric constant p-type
p.eppi = 20*p.epp0;         % Dielectric constant intrinsic
p.eppn = 20*p.epp0;         % Dielectric constant n-type
 
% Energy levels
p.EA = 0;                          % Conduction band energy
p.IP = -1.6;                       % Valence band energy
p.PhiC = -0.15;                    % Cathode workfunction
p.PhiA = -1.45;                    % Anode workfunction
p.Eg = p.EA-p.IP;                      % Band Gap

% Effective density of states and doping concentration and band bending
p.N0 = 1e20;                      % effective Density Of States (eDOS) PVSK (Reenen, 2015) in Cal
p.ND = p.N0*exp((p.PhiC-p.EA)/(p.kB*p.T));  % Based on values for electride workfunctions - ensures self consistency                    
p.NA = p.N0*exp((p.IP-p.PhiA)/(p.kB*p.T));

%%%%% MOBILE ION DEFECT DENSITY %%%%%
p.NI = 1e19;                      % [cm-3]

% Intrinsic Fermi Energy
p.Ei = 0.5*((p.EA+p.IP)+p.kB*p.T*log(p.N0/p.N0));

% Charge Densities
p.ni = p.N0*(exp(-p.Eg/(2*p.kB*p.T)));          % Intrinsic carrier density
p.etln0 = p.N0*exp((p.PhiC-p.EA)/(p.kB*p.T));     % Background density electrons in etl/n-type
p.etlp0 = p.N0*exp((p.IP-p.PhiC)/(p.kB*p.T));     % Background density holes in etl/n-type
p.htln0 = p.N0*exp((p.PhiA-p.EA)/(p.kB*p.T));     % Background density electrons in htl/p-type
p.htlp0 = p.N0*exp((p.IP-p.PhiA)/(p.kB*p.T));     % Background density holes in htl/p-type

% Built in potential based on equilibrium Fermi energies of the conatct
% regions
p.Vbi = (p.kB*p.T/p.q)*log((p.ND*p.NA)/p.ni^2);

% This is a FUDGE at the moment- will only work for symmetric doping
% density... Not sure what the problem is here...

%%%%%%% RECOMBINATION %%%%%%%%%%%
% Radiative recombanation, U = k(np - ni^2)
p.krad = 1e-12;           % [cm3 s-1] Bulk Radiative Recombination coefficient [nominally 1e-10]
p.kradetl = p.krad;         % [cm3 s-1] ETL Radiative Recombination coefficient 
p.kradhtl = p.krad;         % [cm3 s-1] HTL Radiative Recombination coefficient

% SRH recmobination in the contact regions, 
% U = (np-ni^2)/(taun(p+pt) +taup(n+nt))
p.taun_etl = 1e-11;         % [s] SRH time constant for electrons
p.taup_etl = 1e-11;    % [s] SRH time constant for holes
p.taun_htl = 1e-11;        %%%% USE a high value of (e.g.) 1 to switch off
p.taup_htl = 1e-11;    %%%% NOT 0- these variables are in the denominator
p.taun_i = 1e6;
p.taup_i = 1e6;
p.sn = 0;%1e7;            % [cm s-1] electron surface recombination velocity (rate constant for recombination at interface)
p.sp = 0;%sn;             % [cm s-1] hole surface recombination velocity (rate constant for recombination at interface)

% SRH parameters
% se = 1e-15;             % [cm^2] Electron capture cross-section
% v = 1e9;               % [cm/s] Carrier group velocity. Estimates from: http://www.tf.uni-kiel.de/matwis/amat/semi_en/kap_2/advanced/t2_3_1.html
p.Etetl = -0.8;            % ((EA-IP)/2+IP)-0.2;      % Deep trap energy- currently a guess!
p.Ethtl = -0.8;
p.Eti = -0.8;
p.ntetl = p.ni*exp((p.Etetl-p.Ei)/(p.kB*p.T));     % Density of CB electrons when Fermi level at trap state energy
p.ptetl = p.ni*exp((p.Ei-p.Etetl)/(p.kB*p.T));     % Density of VB holes when Fermi level at trap state energy
p.nthtl = p.ni*exp((p.Ethtl-p.Ei)/(p.kB*p.T));     % Density of CB electrons when Fermi level at trap state energy
p.pthtl = p.ni*exp((p.Ei-p.Ethtl)/(p.kB*p.T));     % Density of VB holes when Fermi level at trap state energy
p.nti = p.ni*exp((p.Eti-p.Ei)/(p.kB*p.T));
p.pti = p.ni*exp((p.Ei-p.Eti)/(p.kB*p.T));

% Define geometry of the system (m=0 1D, m=1 cylindrical polar coordinates,
% m=2 spherical polar coordinates).
p.m = 0;

% Define the default relative tolerance for the pdepe solver
% 1e-3 is the default, can be decreased if more precision is needed
p.RelTol = 1e-3;

% Space Charge Region- initial guess required with trial and error better
% than analytical solution
p.wp = 25e-7;  %((-ti*NA*q) + ((NA^0.5)*(q^0.5)*(((ti^2)*NA*q) + (4*eppi*Vbi))^0.5))/(2*NA*q);
p.wn = p.wp;

p.wscr = p.wp + p.ti + p.wn;    % width of space charge region

% Doped region Fermi levels
p.Efnside = - p.Vbi - p.EA + ((p.kB*p.T)/p.q) * log(p.ND/p.N0);
p.Efpside = p.IP - ((p.kB*p.T)/p.q) * log(p.NA/p.N0);

%%%%%%%%%% CURRENT VOLTAGE SCAN %%%%%%%%%%%%
if p.JV == 1
    
    p.tmax = abs(p.Vend- p.Vstart)/p.JVscan_rate;           % Scan time determined by mobility- ensures cell is stable at each point
    p.t0 = 0;
    p.tmesh_type = 1;
    p.tpoints = p.JVscan_pnts;
    
end


% % Pack parameters in to structure 'params'
% varlist = who('*')';
% varstr = strjoin(varlist, ',');
% 
% varcell = who('*')';                    % Store variables names in cell array
% varcell = ['fieldnames', varcell];      % adhere to syntax for v2struct
% 
% params = v2struct(varcell);
params = p; %saves writing params all the time

end
