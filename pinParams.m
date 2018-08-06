function p = pinParams

% Energy levels from Harwell2016a
% Mobilities from 
% Dielectric constants from VanReenen2015a

%Physical constants
p.kB = 8.6173324e-5;    % Boltzmann constant [eV K^-1]
p.T = 300;              % Temperature [K]
p.epp0 = 552434;        % e^2 eV^-1 cm^-1 -Checked (02-11-15)
p.q = 1;                % in e
p.e = 1.61917e-19;      % Charge of an electron in Coulombs for current calculations

p.d = [100e-7, 400e-7, 100e-7];
p.dcum = cumsum(p.d);

% Device Dimensions
p.tp = 100e-7;            % p-type layer thickness
p.pp = 200;               % p-type layer points
p.ti = 400e-7;            % Intrinsic layer thickness
p.pii = 200;              % Intrinsic points called'pii to avoid confusion with the constant 'pi'
p.tn = 100e-7;            % n-type thickness
p.pn = 200;               % n-type points (linear array only)
p.tscr = 30e-7;
p.pscr = 90;
p.tint = 0.2e-7;              % 0.5*Interfacial region thickness Heterojunction
p.pint = 20;               % 0.5*No. of points for interfaciial region
p.te = 10e-7;
p.pepe = 10;
p.tmin = 1e-7;

p.deltax = p.tint/p.pint;        % spacing in the interfacial region

% Parameters for spatial mesh of solution points
p.xmesh_type = 9;     
p.xmax =  p.tp + p.ti + p.tn;   % cm
p.m = 0;  % Define geometry of the system (m=0 1D, m=1 cylindrical polar coordinates,
        % m=2 spherical polar coordinates).

if p.xmesh_type == 1 || p.xmesh_type == 2

    p.x0 = 0;
    
else
    
    p.x0 = p.xmax/1e3;

end

% General Parameters
p.Ana = 1;
p.OC = 0;                 % Toggle open circuit
p.Int = 0;                % Bias Light intensity
p.G0 = 6.25e25 * p.ti;    % Adjusted to give ~ 16mAcm-2
p.tmax = 1e-2;            % Time
p.pulseon = 0;            % Switch pulse on
p.pulselen = 1e-2;        % Transient pulse length
p.pulseint = 0.2;         % Transient pulse intensity- for BL and TM models, 100 mW/cm2 assumed
p.Vapp = 0;               % Applied bias
p.BC = 3;                 % Must be set to one for first solution
p.figson = 1;             % Toggle figures on/off
p.meshx_figon = 0;
p.mesht_figon = 0;
p.side = 1;               % illumination side 1 = EE, 2 = SE
p.calcJ = 0;              % Calculates Currents- slows down solving\
p.JV = 0;

% Solver options
p.RelTol = 1e-3;
p.AbsTol = 1e-6;

% Optical Model 
% 0 = Uniform Generation 
% 1 = Beer-Lamber
% 2 = Transfer Matrix (Stanford)
p.OM  = 1;

p.I0 = 9.5e16;          % Incident photon flux for Beer-Lambert
p.alpha = 6.1e4;        % Effective absorption coefficient

% time parameters
p.tmesh_type = 2;
p.t0 = p.tmax/1e6;
p.tpoints = 400;

%%%%%%%%%% CURRENT VOLTAGE SCAN %%%%%%%%%%%%
% NOT 100% reliable- requires further development and investigation
p.Vstart = 0;
p.Vend = 1.2;
p.JVscan_pnts = 100;

%%%%%%%% General device properties %%%%%%%%%%%%

% Workfunction energies
p.PhiA = 1.9-4.8;
p.PhiC = 1.9-4.4;

%%%%%%%% p-type layer properties %%%%%%%%%%%%%%
p.EA = 1.9 + [-1.9, -3.7, -4.1];
p.IP = 1.9 + [-4.9, -5.3, -7.4];
p.N0 = [1e20, 1e20, 1e20];

% Conversion factors
p.Bn = [(p.N0(1)/p.N0(2))*exp((p.EA(2)-p.EA(1))/(p.kB*p.T)), 1, (p.N0(3)/p.N0(2))*exp((p.EA(2)-p.EA(3))/(p.kB*p.T))];
p.Bp = [(p.N0(1)/p.N0(2))*exp((p.IP(1)-p.IP(2))/(p.kB*p.T)), 1, (p.N0(3)/p.N0(2))*exp((p.IP(3)-p.IP(2))/(p.kB*p.T))];

% For backwards compatibility

% Energy levels (Spiro)
p.EAp = p.EA(1);                      % Conduction band energy (All energies taken from Li et al., 2016, fig 1a)
p.IPp = p.IP(1);                      % Valence band energy
p.Egp = p.EAp - p.IPp;                 % Band Gap

% Energy levels active layer (MAPI)
p.EAi = p.EA(2);                      % Conduction band energy
p.IPi = p.IP(2);                      % Valence band energy
p.Egi = p.EAi - p.IPi;                 % Band Gap

% Energy levels n-type (Ti02)
p.EAn = p.EA(3);                      % Conduction band energy
p.IPn = p.IP(3);                      % Valence band energy
p.Egn = p.EAn - p.IPn;                 % Band Gap

% Effective Density Of States (eDOS) PVSK (Reenen, 2015)
% DIFFERENT EDOS IN DIFFERENT LAYERS UNTESTED!
p.N0p = p.N0(1);                      
p.N0i = p.N0(2); 
p.N0n = p.N0(3);

% Intrinsic carrier concentrations
p.ni = [p.N0p*(exp(-p.Egp/(2*p.kB*p.T))), p.N0i*(exp(-p.Egi/(2*p.kB*p.T))), p.N0n*(exp(-p.Egn/(2*p.kB*p.T)))];

p.nip = p.ni(1);
p.nii = p.ni(2);
p.nin = p.ni(3);

% Intrinsic Fermi Energies
p.Eif = [0.5*((p.EAp+p.IPp)+p.kB*p.T*log(p.N0p/p.N0p)), 0.5*((p.EAi+p.IPi)+p.kB*p.T*log(p.N0i/p.N0i)), 0.5*((p.EAn+p.IPn)+p.kB*p.T*log(p.N0n/p.N0n));];

p.Eifp = p.Eif(1);
p.Eifi = p.Eif(2);
p.Eifn = p.Eif(3);

% Dielectric constants (Van Reenan, 2015)
p.eppp = 3*p.epp0;
p.eppi = 12*p.epp0;
p.eppn = 20*p.epp0;  %20 

% Material Properties
p.mue_p = 20;        % electron mobility (Van Reenan, 2015)
p.muh_p = 20;
p.mue_i = 20;          % electron mobility
p.muh_i = 20;          % hole mobility 
p.mue_n = 20;          % electron mobility
p.muh_n = 20;          % hole mobility   

% Self consistent doping concentrations
p.NA = [p.N0p*exp((p.IP(1)-p.PhiA)/(p.kB*p.T)), 0, 0];
p.ND = [0, 0, p.N0n*exp((p.PhiC-p.EA(3))/(p.kB*p.T))];

% p.NA1 = p.N0p*exp((p.IP(1)-p.PhiA)/(p.kB*p.T));
% p.NA2 = 0;%p.N0n*exp((p.IPn-p.PhiC)/(p.kB*p.T));
% 
% p.ND1 = 0;%p.N0p*exp((p.PhiC-p.EAp)/(p.kB*p.T));
% p.ND2 = p.N0n*exp((p.PhiC-p.EA(3))/(p.kB*p.T));

% Equilibrium Charge Densities
p.n0 = [p.N0(1)*exp((p.PhiA-p.EA(1))/(p.kB*p.T)), p.ni(2), p.N0(3)*exp((p.PhiC-p.EA(3))/(p.kB*p.T))];
p.p0 = [p.N0(1)*exp((p.IP(1)-p.PhiA)/(p.kB*p.T)), p.ni(2), p.N0(3)*exp((p.IP(3)-p.PhiC)/(p.kB*p.T))];

% p.n0p = p.N0p*exp((p.PhiA-p.EA(1))/(p.kB*p.T));     % Background density electrons in htl/p-type
% p.p0p = p.N0p*exp((p.IP(1)-p.PhiA)/(p.kB*p.T));     % Background density holes in htl/p-type
% 
% p.n0n = p.N0n*exp((p.PhiC-p.EA(3))/(p.kB*p.T));     % Background density electrons in etl/n-type
% p.p0n = p.N0n*exp((p.IP(3)-p.PhiC)/(p.kB*p.T));     % Background density holes in etl/n-type

% Doped Equilibrium Fermi Energies
p.E0 = [p.EA(1) + (p.kB*p.T/p.q)*log(p.n0(1)/p.N0(1)), p.Eif(2), p.EA(3) + (p.kB*p.T/p.q)*log(p.n0(3)/p.N0(3))];
% p.Edfp = p.EAp + (p.kB*p.T/p.q)*log(p.n0p/p.N0p);    
% p.Edfn = p.EAn + (p.kB*p.T/p.q)*log(p.n0n/p.N0n); 

% Ion properties
p.NI = 1e19;            % [cm-3] Mobile defect density
p.k_defect_p = 0;      % Defect recombination rate coefficient
p.k_defect_n = 0;
p.mui = 1e-8;

% Vbi based on workfunction difference
p.Vbi = p.PhiC - p.PhiA;

% Vbi for pn junction
% Vbi = (Eifi-Edfp);%(kB*T/q)*log((ND*NA)/nii^2)

% Heterojunction interface
p.DEApi=(p.EAi-p.EAp)/(2*p.tint);
p.DEAin=(p.EAn-p.EAi)/(2*p.tint);
p.DIPpi=(p.IPi-p.IPp)/(2*p.tint);
p.DIPin=(p.IPn-p.IPi)/(2*p.tint);
p.DN0pi=(p.N0i-p.N0p)/(2*p.tint);
p.DN0in=(p.N0n-p.N0i)/(2*p.tint);

% Recombination
p.kradi = 1e-12;      % [cm3 s-1] Bulk Radiative Recombination coefficient [based on SQ limit, 400nm thickness]
p.kradn= 1e-12;       % [cm3 s-1] ETL Radiative Recombination coefficient [based on SQ limit, 100nm thickness]
p.kradp = 1e-12;      % [cm3 s-1] HTL Radiative Recombination coefficient[based on SQ limit, 100nm thickness]

p.taun_i = 1e-8;
p.taup_i = 1e-8;
p.taun_etl = 1e6;       % [s] SRH time constant for electrons
p.taun_htl = 1e6;       % 1e-21;    %%%% USE a high value of (e.g.) 1 to switch off
p.taup_etl = p.taun_etl;    % [s] SRH time constant for holes
p.taup_htl = p.taun_htl;    %%%% NOT 0- these variables are in the denominator

p.sn_r = 1e6;            % [cm s-1] electron surface recombination velocity (rate constant for recombination at interface)
p.sn_l = 1e6;
p.sp_r = 1e6;             % [cm s-1] hole surface recombination velocity (rate constant for recombination at interface)
p.sp_l = 1e6;

% SRH parameters p-type
p.Et_p = p.EAp-p.Egp/2;            %mid-gap trap energy
p.nt_p = p.nip*exp((p.Et_p-p.Eifp)/(p.kB*p.T));     % Density of CB electrons when Fermi level at trap state energy
p.pt_p = p.nip*exp((p.Eifp-p.Et_p)/(p.kB*p.T));     % Density of VB holes when Fermi level at trap state energy

% SRH parameters intrinsic
p.Et_i = p.EAi-p.Egi/2;            %mid-gap trap energy
p.nt_i = p.nii*exp((p.Et_i-p.Eifi)/(p.kB*p.T));     % Density of CB electrons when Fermi level at trap state energy
p.pt_i = p.nii*exp((p.Eifi-p.Et_i)/(p.kB*p.T));     % Density of VB holes when Fermi level at trap state energy

% SRH parameters p-type
p.Et_n = p.EAn-p.Egn/2;            %mid-gap trap energy
p.nt_n = p.nin*exp((p.Et_n-p.Eifn)/(p.kB*p.T));     % Density of CB electrons when Fermi level at trap state energy
p.pt_n = p.nin*exp((p.Eifn-p.Et_n)/(p.kB*p.T));     % Density of VB holes when Fermi level at trap state energy

% Space Charge Region - better to guess than use analytical solutions
p.wp = 10e-7;  %((-ti*NA*q) + ((NA^0.5)*(q^0.5)*(((ti^2)*NA*q) + (4*eppi*Vbi))^0.5))/(2*NA*q);
p.wn = 10e-7;

p.wscr = p.wp + p.ti + p.wn;

end


%%
%{
%%%%%%%% TRANSIENT SETTINGS %%%%%%%%
if p.pulseon == 1 && p.OC == 1        % Record length for TPV
    
    p.pulselen = 60e-6;       % Transient pulse length
    p.pulseint = 2;          % Transient pulse intensity- for BL and TM models, 100 mW/cm2 assumed
    p.pulsestart = 1e-7;     % Time before pulse- required for baseline
    p.tmesh_type = 3;
    p.tmax = 100e-6;
    p.t0 = 0;
    p.tpoints = 400;
    p.deltat = p.tmax/(1e4*p.tpoints);

end    

if p.pulseon == 1 && p.OC == 0        % Record length for TPC
    
    p.pulselen = 10e-6;       % Transient pulse length
    p.pulseint = 0.2;            % Transient pulse intensity- for BL and TM models, 100 mW/cm2 assumed
    p.pulsestart = 1e-6;
    p.tmesh_type = 1;
    p.tmax = 40e-6;  
    p.t0 = 0;%tmax/1000;
    p.tpoints = 400;
    p.deltat = p.tmax/(1e4*p.tpoints);
    p.calcJ = 1;
    
end
%}