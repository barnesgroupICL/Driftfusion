% Class defining the properties and methods for the simulation parameters
classdef pc
    % Parameters Class
    
    properties (Constant)
        % Set device properties
        % These cannot be altered for subsequent solutions once set
        
        % Physical constants
        
        kB = 8.617330350e-5;       % Boltzmann constant [eV K^-1]
        epp0 = 552434;           % Epsilon_0 [e^2 eV^-1 cm^-1] - Checked (02-11-15)
        q = 1;                   % Charge of the species in units of e.
        e = 1.61917e-19;         % Elementary charge in Coulombs.
        
        % Device Dimensions [cm]
        d = [22e-7, 510e-7, 60e-7];       % Layer thickness array
        parr = [10, 200, 40];             % Spatial mesh points array
        
        dcum = cumsum(pc.d);                % cumulative thickness        
                 
        dint = 2e-7;        % 0.5x Interfacial region thickness (x_mesh_type = 3), this is related to Space Charge Region, read below for wp and wn parameters
        pint = 20;          % 0.5x Interfacial points (x_mesh_type = 3)
        dscr = 20e-7;       % Approx space charge region thickness
        pscr = 40;          % No of points in space charge region
        pepe = 20;          % electrode interface points
        te = 10e-7;         % electrode interface thickness
        dmin = 1e-7;        % start value for log meshes
        deltax = 1e-7;
        
        % Define spatial cordinate system
        % m=0 cartesian
        % m=1 cylindrical polar coordinates
        % m=2 spherical polar coordinates
        m = 0;
        
    end
    
    properties
        % Sets the default values
        
        % Temperature [K]
        T = 300;
        
        %%%%%%% GENERAL CONTROL PARAMETERS %%%%%%%%%%
        OC = 0;                 % Closed circuit = 0, Open Circuit = 1
        Int = 0;                % Bias Light intensity (Suns Eq.)
        G0 = 2.6409e+21;        % Uniform generation rate @ 1 Sun
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
        OM = 0;
        
        %%%%%%% MESHES %%%%%%
        % xmesh_type specification - see xmesh_gen
        xmesh_type = 3;
        
        % Parameters for time mesh
        tmax = 1e-12;               % Max time value
        t0 = 1e-16;                 % Initial log mesh time value
        tpoints = 100;              % Number of time points
        tmesh_type = 2;             % Mesh type- for use with meshgen_t
        
        %%%%%%%%%%% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%
        %% Layer description- currently for optical properties only
        % See Index of Refraction library for choices- names must be exactly the same but '_n', '_k' should be omitted
        stack = {'PS', 'MAPICl', 'PCBM'}
        
        %% Energy levels
        %PEDOT:PSS        
        EA = [-3.5, -3.8, -3.8];%   %1.9 + [-1.9, -3.7, -4.1];
        IP = [-5.1, -5.4, -6.2];%    %1.9 + [-4.9, -5.3, -7.4]; 
        
        % PCBM: Sigma Aldrich https://www.sigmaaldrich.com/technical-documents/articles/materials-science/organic-electronics/pcbm-n-type-semiconductors.html
        
        %% Equilibrium Fermi energies - defines doping density
        E0 = [-4.3, -4.6, -3.9];
        %E0 = [-5.25, -4.6, -3.95];
        
        % Workfunction energies
        PhiA = -5.0;    %-1.4;    %
        PhiC = -3.9;    %-0.6;    %
        
        % Effective Density Of States
        % DIFFERENT eDOS IN DIFFERENT LAYERS AS YET UNTESTED!
        N0 = [1e19, 6e18, 1e19];
        % PEDOT eDOS: https://aip.scitation.org/doi/10.1063/1.4824104
        % MAPI eDOS: F. Brivio, K. T. Butler, A. Walsh and M. van Schilfgaarde, Phys. Rev. B, 2014, 89, 155204.
        % PCBM eDOS:
                
        %%%%% MOBILE ION DEFECT DENSITY %%%%%
        NI = 1e18;                      % [cm-3] A. Walsh, D. O. Scanlon, S. Chen, X. G. Gong and S.-H. Wei, Angewandte Chemie, 2015, 127, 1811.
        a_max = 1.21e22;                % [cm-3] max density of iodide sites- P. Calado thesis
        
        % Mobilities
        mue = [1e-4, 20, 1e-3];         % electron mobility [cm2V-1s-1]
        muh = [1e-4, 20, 1e-3];         % hole mobility [cm2V-1s-1]
        mui = 1e-10;                   % ion mobility [cm2V-1s-1]
        % PTPD h+ mobility: https://pubs.rsc.org/en/content/articlehtml/2014/ra/c4ra05564k
        % PEDOT e- mobility: 0.01 cm2V-1s-1 https://aip.scitation.org/doi/10.1063/1.4824104
        
        % Dielectric constants
        epp = [4,23,4];
        
        %%%%%%% RECOMBINATION %%%%%%%%%%%
        % Radiative recombination, U = k(np - ni^2)
        krad = [8.3e-11, 1.0e-11, 1.0e-11, 1.0e-11, 6.8e-11];  % [cm3 s-1] Radiative Recombination coefficient
        % p-type/pi-interface/intrinsic/in-interface/n-type
        
        taun = [1e-9, 1e-6, 1e-6];           % [s] SRH time constant for electrons
        taup = [1e-9, 1e-6, 1e-6];           % [s] SRH time constant for holes
        
        % Surface recombination and extraction coefficients
        sn_r = 1e8;            % [cm s-1] electron surface recombination velocity (rate constant for recombination at interface)
        sn_l = 1e8;
        sp_r = 1e8;             % [cm s-1] hole surface recombination velocity (rate constant for recombination at interface)
        sp_l = 1e8;
        
        % Defect recombination rate coefficient - not currently used
        k_defect_p = 0;
        k_defect_n = 0;
        
        %% Pulse settings
        laserlambda = 638;      % Pulse wavelength (Currently implemented with OM2 (Transfer Matrix only)
        pulselen = 1e-6;        % Transient pulse length
        pulsepow = 10;          % Pulse power [mW cm-2] OM2 (Transfer Matrix only)
        pulsestart = 1e-7;
        
        %% Current voltage scan parameters
        Vstart = 0;             % Initial scan point
        Vend = 1.2;             % Final scan point
        JVscan_rate = 1;        % JV scan rate (Vs-1)
        JVscan_pnts = 100;
        
        % xmax changes dependent on SC or OC - this may be changed in
        % future versions
        xmax = sum(pc.d);
        
        %% Dynamically created variables
        genspace = [];
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
        
        dEAdx
        dIPdx
        dN0dx
        Bn
        Bp
        Eg
        Eif
        Et
        NA
        ND
        Vbi
        n0
        nleft
        nright
        ni
        nt          % Density of CB electrons when Fermi level at trap state energy
        p0
        pleft
        pright
        pt           % Density of VB holes when Fermi level at trap state energy
        wn
        wp
        wscr            % Space charge region width
        x0              % Initial spatial mesh value
        
    end
    
    
    methods
        
        function params = pc
            % paramsStruct Constructor for paramsStruct.
            %
            % Warn if tmesh_type is not correct
            
            if ~ any([1 2 3 4] == params.tmesh_type)
                warning('PARAMS.tmesh_type should be an integer from 1 to 4 inclusive. MESHGEN_T cannot generate a mesh if this is not the case.')
            end
            
            % Warn if xmesh_type is not correct
            if ~ any(1:1:10 == params.xmesh_type)
                warning('PARAMS.xmesh_type should be an integer from 1 to 9 inclusive. MESHGEN_X cannot generate a mesh if this is not the case.')
            end
            
            for i = 1:length(params.ND)
                if params.ND(i) >= params.N0(i) || params.NA(i) >= params.N0(i)
                    msg = 'Doping density must be less than eDOS. For consistent values ensure electrode workfunctions are within the band gap and check expressions for doping density in Dependent variables.';
                    error(msg);
                end
            end
        end
        
        function params = set.xmesh_type(params, value)
            % SET.xmesh_type Check if xmesh_type is an integer from 1 to 4.
            %
            %   SET.xmesh_type(PARAMS, VALUE) checks if VALUE is an integer
            %   from 1 to 4, and if so, changes PARAMS.xmesh_type to VALUE.
            %   Otherwise, a warning is shown. Runs automatically whenever
            %   xmesh_type is changed.
            if any(1:1:9 == value)
                params.xmesh_type = value;
            else
                error('PARAMS.xmesh_type should be an integer from 1 to 9 inclusive. MESHGEN_X cannot generate a mesh if this is not the case.')
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
        
        function params = set.ND(params, value)
            for i = 1:length(params.ND)
                if value(i) >= params.N0(i)
                    error('Doping density must be less than eDOS. For consistent values ensure electrode workfunctions are within the band gap.')
                end
            end
        end
        
        function params = set.NA(params, value)
            for i = 1:length(params.ND)
                if value(i) >= params.N0(i)
                    error('Doping density must be less than eDOS. For consistent values ensure electrode workfunctions are within the band gap.')
                end
            end
        end
               
        % Band gap energies
        function value = get.Eg(params)
            
            value = params.EA - params.IP;
            
        end
        
        % Vbi based on difference in boundary workfunctions
        function value = get.Vbi(params)
            
            value = params.PhiC - params.PhiA;
            
        end
        
        % Boltzmann conversion factor electrons
        function value = get.Bn(params)
                
%             value = [(params.N0(1)/params.N0(2))*exp((params.EA(2)-params.EA(1))/(params.kB*params.T)), 1,...
%                 (params.N0(3)/params.N0(2))*exp((params.EA(2)-params.EA(3))/(params.kB*params.T))];
               value = [1,1,1];
        end
        
        % Boltzmann conversion factor holes
        function value = get.Bp(params)
            
            %Bp = [1, (params.N0(2)/params.N0(1))*exp((params.IP(2)-params.IP(1))/(params.kB*params.T)), (params.N0(3)/params.N0(1))*exp((params.IP(3)-params.IP(1))/(params.kB*params.T))];
            %Bp = [(params.N0(1)/params.N0(3))*exp((params.IP(1)-params.IP(3))/(params.kB*params.T)), (params.N0(2)/params.N0(3))*exp((params.IP(2)-params.IP(3))/(params.kB*params.T)), 1];
%             value = [(params.N0(1)/params.N0(2))*exp((params.IP(1)-params.IP(2))/(params.kB*params.T)), 1,...
%                 (params.N0(3)/params.N0(2))*exp((params.IP(3)-params.IP(2))/(params.kB*params.T))];
%             
                value = [1,1,1];
        end
        
        % Intrinsic Fermi Energies
        function value = get.Eif(params)
            
            value = [0.5*((params.EA(1)+params.IP(1))+params.kB*params.T*log(params.N0(1)/params.N0(1))),...
                0.5*((params.EA(2)+params.IP(2))+params.kB*params.T*log(params.N0(2)/params.N0(2))),...
                0.5*((params.EA(3)+params.IP(3))+params.kB*params.T*log(params.N0(3)/params.N0(3)))];
            
        end
        
%         % Doped Equilibrium Fermi Energies
%         function value = get.E0(params)
%             
%             value = [params.EA(1) + (params.kB*params.T/params.q)*log(params.n0(1)/params.N0(1)),...
%                 params.EA(2) + (params.kB*params.T/params.q)*log(params.n0(2)/params.N0(2)),...
%                 params.EA(3) + (params.kB*params.T/params.q)*log(params.n0(3)/params.N0(3))];
%             
%         end
        % Conduction band gradients at interfaces
        function value = get.dEAdx(params)
            
            value = [(params.EA(2)-params.EA(1))/(2*params.dint), (params.EA(3)-params.EA(2))/(2*params.dint)];
         
        end
                
        %Valence band gradients at interfaces
        function value = get.dIPdx(params)
            
            value = [(params.IP(2)-params.IP(1))/(2*params.dint), (params.IP(3)-params.IP(2))/(2*params.dint)];
         
        end
        
        % eDOS gradients at interfaces 
        function value = get.dN0dx(params)
            
            value = [(params.N0(2)-params.N0(1))/(2*params.dint), (params.N0(3)-params.N0(2))/(2*params.dint)];
         
        end

        % Donor densities
        function value = get.ND(params)
            
            value = [0, 0, params.N0(3)*exp((params.E0(3)-params.EA(3))/(params.kB*params.T))];
            
        end
        % Donor densities
        function value = get.NA(params)
            
            value = [params.N0(1)*exp((params.IP(1)-params.E0(1))/(params.kB*params.T)), 0, 0];
            
        end
        
        % Intrinsic carrier densities
        function value = get.ni(params)
            
            value = [params.N0(1)*(exp(-params.Eg(1)/(2*params.kB*params.T))),...
                params.N0(2)*(exp(-params.Eg(2)/(2*params.kB*params.T))),...
                params.N0(3)*(exp(-params.Eg(3)/(2*params.kB*params.T)))];
            
        end
        % Equilibrium electron diensities
        function value = get.n0(params)
            
            value = [params.N0(1)*exp((params.E0(1)-params.EA(1))/(params.kB*params.T)),...
                params.N0(2)*exp((params.E0(2)-params.EA(2))/(params.kB*params.T)),...
                params.N0(3)*exp((params.E0(3)-params.EA(3))/(params.kB*params.T))];     % Background density electrons in htl/p-type
            
        end
              
        
        function value = get.p0(params)
            
            value = [params.N0(1)*exp((params.IP(1)-params.E0(1))/(params.kB*params.T)),...
                params.N0(2)*exp((params.IP(2)-params.E0(2))/(params.kB*params.T)),...
                params.N0(3)*exp((params.IP(3)-params.E0(3))/(params.kB*params.T))];     % Background density holes in htl/p-type
            
        end
        
        
        % Boundary electron and hole densities - based on the workfunction
        % of the electrodes
        function value = get.nleft(params)
            
            value = params.N0(1)*exp((params.PhiA-params.EA(1))/(params.kB*params.T));
        
        end
        
                % Boundary electron and hole densities
        function value = get.nright(params)
            
            value = params.N0(3)*exp((params.PhiC-params.EA(3))/(params.kB*params.T));
        
        end
        
                % Boundary electron and hole densities
        function value = get.pleft(params)
            
            value = params.N0(1)*exp((params.IP(1)-params.PhiA)/(params.kB*params.T));
        
        end
        
        function value = get.pright(params)
            
            value = params.N0(3)*exp((params.IP(3)-params.PhiC)/(params.kB*params.T));
        
        end
        
        %% SRH parameters
        
        % SRH trap energies- currently set to mid gap - always set
        % respective to energy levels to avoid conflicts
        % U = (np-ni^2)/(taun(p+pt) +taup(n+nt))
        function value = get.Et(params)
            
            value = [params.EA(2)-(params.Eg(2)/2), params.EA(2)-(params.Eg(2)/2), params.EA(2)-(params.Eg(2)/2)];
            
        end
        
        % nt_p - Density of CB electrons when Fermi level at trap state energy
        function value = get.nt(params)
            
            value = params.ni.*exp((params.Et-params.Eif)/(params.kB*params.T));       % cm
            
        end
        
        % pt_p - Density of VB holes when Fermi level at trap state energy
        function value = get.pt(params)
            
            value =  params.ni.*exp((params.Eif-params.Et)/(params.kB*params.T));      % cm
            
        end
        
        %% Space charge layer widths
        function value = get.wp (params)
            
            value = ((-params.d(2)*params.NA(1)*params.q) + ((params.NA(1)^0.5)*(params.q^0.5)*(((params.d(2)^2)*params.NA(1)*params.q) + (4*params.epp(2)*params.Vbi))^0.5))/(2*params.NA(1)*params.q);
            
        end
        
        function value = get.wn (params)
            
            value = ((-params.d(2)*params.ND(3)*params.q) + ((params.ND(3)^0.5)*(params.q^0.5)*(((params.d(2)^2)*params.ND(3)*params.q) + (4*params.epp(2)*params.Vbi))^0.5))/(2*params.ND(3)*params.q);
            
        end
        
        % wscr - space charge region width
        function value = get.wscr(params)
            
            value = params.wp + params.d(2) + params.wn;         % cm
            
        end
        
        %%
        
    end
end
