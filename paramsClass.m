classdef paramsClass
    
    properties (Constant)
        % Set device properties
        % These cannot be altered for subsequent solutions once set
        
        % Physical constants
        
        kB = 8.6173324e-5;       % Boltzmann constant [eV K^-1]
        T = 300;                 % Temperature [K]
        epp0 = 552434;           % Epsilon_0 [e^2 eV^-1 cm^-1] - Checked (02-11-15)
        q = 1;                   % Charge of the species in units of e.
        e = 1.61917e-19;         % Elementary charge in Coulombs.
        
        % Device Dimensions [cm]
        d = [100e-7, 400e-7, 100e-7];       % Layer dimensions array
        dmax = sum(paramsClass.d);        
        
        tp = 200e-7;        % p-type layer thickness
        pp = 100;            % p-type layer points
        ti = 400e-7;        % Intrinsic layer thickness; this affects the wp parameter
        pii = 200;          % Intrinsic points
        tn = 200e-7;        % n-type thickness
        pn = 100;            % n-type points
        tint = 2e-7;        % 0.5x Interfacial region thickness (x_mesh_type = 3), this is related to Space Charge Region, read below for wp and wn parameters
        pint = 40;          % 0.5x Interfacial points (x_mesh_type = 3)
        tscr = 50e-7;       % Approx space charge region thickness
        pscr = 50;          % No of points in space charge region
        pepe = 20;          % electrode interface points
        te = 10e-7;         % electrode interface thickness  
        tmin = 1e-7;        % start value for log meshes
        
        % Define spatial cordinate system
        % m=0 cartesian
        % m=1 cylindrical polar coordinates
        % m=2 spherical polar coordinates
        m = 0;
        
        %% Energy levels
        EA = [-1.9, -3.7, -4.1];
        IP = [-4.9, -5.3, -7.4];
        Eg = [paEA(1) - IP
        
        % Workfunction energies
        PhiA = -4.8;    % -5.1; Gold (Van Reenan, 2015)
        PhiC = -4.4;    % FTO (Van Reenan, 2015)
        
        % For backwards compatibility
        % Energy levels (HTL)
        EAp = paramsClass.EA(1);                      % Conduction band energy (All energies taken from Li et al., 2016, fig 1a)
        IPp = paramsClass.IP(1);                      % Valence band energy
        Egp = paramsClass.EA(1) - paramsClass.IP(1);                 % Band Gap

        % Energy levels active layer (Intrinsic)
        EAi = paramsClass.EA(2);                      % Conduction band energy
        IPi = paramsClass.IP(2);                      % Valence band energy
        Egi = paramsClass.EA(2) - paramsClass.IP(2);                 % Band Gap

        % Energy levels n-type (ETL)
        EAn = paramsClass.EA(3);                      % Conduction band energy
        IPn = paramsClass.EA(3);                      % Valence band energy
        Egn = paramsClass.EA(3) - paramsClass.IP(3);                 % Band Gap
        
        % Vbi based on workfunction difference
        Vbi = paramsClass.PhiC - paramsClass.PhiA;
              
        % SRH trap energies- currently set to mid gap - always set
        % respective to energy levels to avoid conflicts
        % U = (np-ni^2)/(taun(p+pt) +taup(n+nt))
        Et_p = paramsClass.EAp-(paramsClass.Egp/2);           % SRH trap energy p-type
        Et_i = paramsClass.EAi-(paramsClass.Egi/2);           % SRH trap energy intrinsic
        Et_n = paramsClass.EAn-(paramsClass.Egn/2);           % SRH trap energy n-type
                
        % Effective Density Of States (eDOS) PVSK (Reenen, 2015)
        % DIFFERENT eDOS IN DIFFERENT LAYERS UNTESTED!
        N0 = [1e20, 1e20, 1e20];
        
        % For backwards compatibility
        N0p = paramsClass.N0(1);
        N0i = paramsClass.N0(2);
        N0n = paramsClass.N0(3); 
                
         % Conversion factors
        Bn = [(paramsClass.N0(1)/paramsClass.N0(2))*exp((paramsClass.EA(2)-paramsClass.EA(1))/(paramsClass.kB*paramsClass.T)), 1, (paramsClass.N0(3)/paramsClass.N0(2))*exp((paramsClass.EA(2)-paramsClass.EA(3))/(paramsClass.kB*paramsClass.T))];
        Bp = [(paramsClass.N0(1)/paramsClass.N0(2))*exp((paramsClass.IP(1)-paramsClass.IP(2))/(paramsClass.kB*paramsClass.T)), 1, (paramsClass.N0(3)/paramsClass.N0(2))*exp((paramsClass.IP(3)-paramsClass.IP(2))/(paramsClass.kB*paramsClass.T))];
        
        % Intrinsic Fermi Energies
        Eifp = 0.5*((paramsClass.EAp+paramsClass.IPp)+paramsClass.kB*paramsClass.T*log(paramsClass.N0p/paramsClass.N0p));
        Eifi = 0.5*((paramsClass.EAi+paramsClass.IPi)+paramsClass.kB*paramsClass.T*log(paramsClass.N0i/paramsClass.N0i));
        Eifn = 0.5*((paramsClass.EAn+paramsClass.IPn)+paramsClass.kB*paramsClass.T*log(paramsClass.N0n/paramsClass.N0n));
        
        % Doped Equilibrium Fermi Energies
        Edfp = paramsClass.EAp + (paramsClass.kB*paramsClass.T/paramsClass.q)*log(paramsClass.n0p/paramsClass.N0p);    
        Edfn = paramsClass.EAn + (paramsClass.kB*paramsClass.T/paramsClass.q)*log(paramsClass.n0n/paramsClass.N0n); 
        
        % Self consistent doping concentrations - HTL and ETL in
        % equilibrium with electrodes- does not need to be the case
        %NA = paramsClass.N0p*exp((paramsClass.IPp-paramsClass.Edfp)/(paramsClass.kB*paramsClass.T));
        %ND = paramsClass.N0n*exp((paramsClass.Edfn-paramsClass.EAn)/(paramsClass.kB*paramsClass.T));
        
        % Boundary Charge Densities
        n0p = paramsClass.N0(1)*exp((paramsClass.PhiA-paramsClass.EA(1))/(paramsClass.kB*paramsClass.T));     % Background density electrons in htl/p-type
        p0p = paramsClass.N0(1)*exp((paramsClass.IP(1)-paramsClass.PhiA)/(paramsClass.kB*paramsClass.T));     % Background density holes in htl/p-type

        n0n = paramsClass.N0(3)*exp((paramsClass.PhiC-paramsClass.EA(3))/(paramsClass.kB*paramsClass.T));     % Background density electrons in etl/n-type
        p0n = paramsClass.N0(3)*exp((paramsClass.IP(3)-paramsClass.PhiC)/(paramsClass.kB*paramsClass.T));     % Background density holes in etl/n-type
      
        % Intrinsic carrier densities
        nip = paramsClass.N0p*(exp(-paramsClass.Egp/(2*paramsClass.kB*paramsClass.T)));     % Intrinsic carrier density p-type
        nii = paramsClass.N0i*(exp(-paramsClass.Egi/(2*paramsClass.kB*paramsClass.T)));     % Intrinsic carrier density
        nin = paramsClass.N0n*(exp(-paramsClass.Egn/(2*paramsClass.kB*paramsClass.T)));     % Intrinsic carrier density p-type

        %%%%% MOBILE ION DEFECT DENSITY %%%%%
        NI = 1e19;                      % [cm-3]
        
        % Self consistent doping concentrations
        NA1 = paramsClass.N0p*exp((paramsClass.IP(1)-paramsClass.PhiA)/(paramsClass.kB*paramsClass.T));
        NA2 = 0;%paramsClass.N0n*exp((paramsClass.IPn-paramsClass.PhiC)/(paramsClass.kB*paramsClass.T));

        ND1 = 0;%paramsClass.N0p*exp((paramsClass.PhiC-paramsClass.EAp)/(paramsClass.kB*paramsClass.T));
        ND2 = paramsClass.N0n*exp((paramsClass.PhiC-paramsClass.EA(3))/(paramsClass.kB*paramsClass.T));
        
    end
    
    
    properties
    % Sets the default values
       
        %%%%%%% GENERAL CONTROL PARAMETERS %%%%%%%%%%
        OC = 0;                 % Closed circuit = 0, Open Circuit = 1
        Int = 0;                % Bias Light intensity (Suns Eq.)
        G0 = 2.5e21;            % Uniform generation rate @ 1 Sun
        pulseon = 0;            % Switch pulse on TPC or TPV
        Vapp = 0;               % Applied bias
        BC = 3;                 % Boundary Conditions. Must be set to one for first solution
        figson = 1;             % Toggle figures on/off
        meshx_figon = 0;        % Toggles x-mesh figures on/off
        mesht_figon = 0;        % Toggles t-mesh figures on/off
        side = 1;               % illumination side 1 = EE, 2 = SE
        calcJ = 0;              % Calculates Currents- slows down solving calcJ = 1, calculates DD currents at every position
        mobset = 1;             % Switch on/off electron hole mobility- MUST BE SET TO ZERO FOR INITIAL SOLUTION
        JV = 0;                 % Toggle run JV scan on/off
        Ana = 1;                % Toggle on/off analysis
        
        % OM = Optical Model
        % Current only uniform generation functionality is avaiable- this will be
        % updated in future versions.
        % 0 = Uniform Generation
        OM  = 0;
        
        %%%%%%% MESHES %%%%%%
        % xmesh_type specification
        xmesh_type = 1;
        
        % Parameters for time mesh
        tmax = 1e-12;               % Max time value
        t0 = 1e-16;                 % Initial log mesh time value
        tpoints = 200;              % Number of time points
        tmesh_type = 2;             % Mesh type- for use with meshgen_t
        
        %%%%%%%%%%% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%
        
        % Mobilities
        mue_p = 20;         % electron mobility p-type
        muh_p = 20;         % hole mobility p-type
        mue_i = 20;         % electron mobility intrinsic
        muh_i = 20;         % hole mobility intrinsic
        mue_n = 20;         % electron mobility n-type
        muh_n = 20;         % hole mobility n-type
        
        mui = 1e-10;        % ion mobility
        
        % Dielectric constants
        eppp = 20*paramsClass.epp0;         % Dielectric constant p-type
        eppi = 20*paramsClass.epp0;         % Dielectric constant intrinsic
        eppn = 20*paramsClass.epp0;         % Dielectric constant n-type
        
        %%%%%%% RECOMBINATION %%%%%%%%%%%
        % Radiative recombination, U = k(np - ni^2)
        kradi = 1e-12;           % [cm3 s-1] Bulk Radiative Recombination coefficient [nominally 1e-10]
        kradp = 1e-12;         % [cm3 s-1] ETL Radiative Recombination coefficient
        kradn = 1e-12;         % [cm3 s-1] HTL Radiative Recombination coefficient
                       
        taun_etl = 1e6;           % [s] SRH time constant for electrons
        taup_etl = 1e6;           % [s] SRH time constant for holes
        taun_htl = 1e6;           %%%% USE a high value of (e.g.) 1 to switch off
        taup_htl = 1e6;           %%%% NOT 0- these variables are in the denominator
        taun_i = 1e-8;
        taup_i = 1e-8;
        
        % Surface recombination and extraction coefficients
        sn_r = 1e8;            % [cm s-1] electron surface recombination velocity (rate constant for recombination at interface)
        sn_l = 1e8;
        sp_r = 1e8;             % [cm s-1] hole surface recombination velocity (rate constant for recombination at interface)
        sp_l = 1e8;

        % Defect recombination rate coefficient - not currently used
        k_defect_p = 0;      
        k_defect_n = 0;
        
        %% Pulse settings
        pulselen = 1e-6;            % Transient pulse length
        pulseint = 0.2;             % Transient pulse intensity- for BL and TM models, 100 mW/cm2 assumed
        pulsestart = 1e-7;
        
        %% Current voltage scan parameters
        Vstart = 0;             % Initial scan point
        Vend = 1.2;             % Final scan point
        JVscan_rate = 1;        % JV scan rate (Vs-1)
        JVscan_pnts = 600;
        
        % xmax changes dependent on SC or OC - this may be changed in
        % future versions
        xmax = paramsClass.tp + paramsClass.ti + paramsClass.tn;
        
        %% Dynamically created variables
        
        x = [];
        t = [];
        xpoints = [];
        
        % Define the default relative tolerance for the pdepe solver
        % 1e-3 is the default, can be decreased if more precision is needed
        % Solver options
        RelTol = 1e-3;
        AbsTol = 1e-6;
        
    end
    
    
    %%  Properties whose values depend on other properties (see 'get' methods).
    properties (Dependent)
    
        wn
        wp
        wscr            % Space charge region width
        x0              % Initial spatial mesh value
        nt_p           % Density of CB electrons when Fermi level at trap state energy
        pt_p           % Density of VB holes when Fermi level at trap state energy
        nt_i 
        pt_i
        nt_n          % Density of CB electrons when Fermi level at trap state energy
        pt_n           % Density of VB holes when Fermi level at trap state energy
        
    end
    
        
    methods
        
        function params = paramsClass
            % paramsStruct Constructor for paramsStruct.
            %
            % Warn if tmesh_type is not correct
            
            if ~ any([1 2 3 4] == params.tmesh_type)
                warning('PARAMS.tmesh_type should be an integer from 1 to 4 inclusive. MESHGEN_T cannot generate a mesh if this is not the case.')
            end
            
            % Warn if xmesh_type is not correct
            if ~ any([1 2 3 4] == params.xmesh_type)
                warning('PARAMS.xmesh_type should be an integer from 1 to 4 inclusive. MESHGEN_X cannot generate a mesh if this is not the case.')
            end
            
        end
        
        function params = set.xmesh_type(params, value)
            % SET.xmesh_type Check if xmesh_type is an integer from 1 to 4.
            %   
            %   SET.xmesh_type(PARAMS, VALUE) checks if VALUE is an integer
            %   from 1 to 4, and if so, changes PARAMS.xmesh_type to VALUE.
            %   Otherwise, a warning is shown. Runs automatically whenever
            %   xmesh_type is changed.
            if any([1 2 3 4] == value)
                params.xmesh_type = value;
            else
                error('PARAMS.xmesh_type should be an integer from 1 to 4 inclusive. MESHGEN_X cannot generate a mesh if this is not the case.')
            end
        end

        function params = set.tmesh_type(params, value)
            % SET.tmesh_type Check if xmesh_type is an integer from 1 to 2.
            % 
            %   SET.tmesh_type(PARAMS, VALUE) checks if VALUE is an integer
            %   from 1 to 2, and if so, changes PARAMS.tmesh_type to VALUE.
            %   Otherwise, a warning is shown. Runs automatically whenever
            %   tmesh_type is changed.
            if any([1 2 3 4] == value)
                params.tmesh_type = value;
            else
                error('PARAMS.tmesh_type should be an integer from 1 to 4 inclusive. MESHGEN_T cannot generate a mesh if this is not the case.')
            end
        end
        
       % set linear mesh for time if JV =1 
       function value = get.tmesh_type(params)      
            
           if params.JV == 1
               
               value = 1;
                
           else
               
               value = params.tmesh_type;
                
           end
           
       end
       
      
%        % xmax - space charge region width
%        function value = get.xmax(params)      
%             
%             value = params.tp + params.ti + params.tn;       % cm
%        
%        end
       
       %% SRH parameters       
       % nt_p - Density of CB electrons when Fermi level at trap state energy
       function value = get.nt_p(params)      
            
            value = params.nip*exp((params.Et_p-params.Eifp)/(params.kB*params.T));       % cm
       
       end
       
       % pt_p - Density of VB holes when Fermi level at trap state energy
       function value = get.pt_p(params)      
            
            value =  params.nip*exp((params.Eifp-params.Et_p)/(params.kB*params.T));      % cm
       
       end
       
       % nt_i
       function value = get.nt_i(params)      
            
            value =  params.nii*exp((params.Et_i-params.Eifi)/(params.kB*params.T));    % cm
       
       end
       
       % pt_i
       function value = get.pt_i(params)
            
            value = params.nii*exp((params.Eifi-params.Et_i)/(params.kB*params.T));     % cm
       
       end        
       
       % nt_n - Density of CB electrons when Fermi level at trap state energy
       function value = get.nt_n (params)
            
            value = params.nin*exp((params.Et_n-params.Eifn)/(params.kB*params.T));         % cm
       
       end        
                 
       % nt_n - Density of VB holes when Fermi level at trap state energy
       function value = get.pt_n (params)
            
            value = params.nin*exp((params.Eifn-params.Et_n)/(params.kB*params.T));      % cm
       
       end
       
       %% Space charge layer widths
       function value = get.wp (params)
           
           value = ((-params.ti*params.NA*params.q) + ((params.NA^0.5)*(params.q^0.5)*(((params.ti^2)*params.NA*params.q) + (4*params.eppi*params.Vbi))^0.5))/(2*params.NA*params.q);
           
       end
       
       function value = get.wn (params)
       
           value = ((-params.ti*params.NA*params.q) + ((params.NA^0.5)*(params.q^0.5)*(((params.ti^2)*params.NA*params.q) + (4*params.eppi*params.Vbi))^0.5))/(2*params.NA*params.q);
       
       end
       
              % wscr - space charge region width
       function value = get.wscr(params)      
            
            value = params.wp + params.ti + params.wn;         % cm
       
       end
       
       %%
           
    end
end
