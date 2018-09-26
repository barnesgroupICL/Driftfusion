% Class defining the properties and methods for the simulation parameters
classdef pc
    % Parameters Class
    
    properties (Constant)
        
        % Physical constants
        kB = 8.617330350e-5;     % Boltzmann constant [eV K^-1]
        epp0 = 552434;           % Epsilon_0 [e^2 eV^-1 cm^-1] - Checked (02-11-15)
        q = 1;                   % Charge of the species in units of e.
        e = 1.61917e-19;         % Elementary charge in Coulombs.
       
    end
    
    properties
        
        % Temperature [K]
        T = 300;
        
        % Device Dimensions [cm]
        d = [100e-7, 400e-7, 100e-7];       % Layer thickness array
        parr = [50, 100, 50];             % Spatial mesh points array
                         
        dint = 2e-7;        % 0.5x Interfacial region thickness (x_mesh_type = 3), this is related to Space Charge Region, read below for wp and wn parameters
        pint = 10;          % 0.5x Interfacial points (x_mesh_type = 3)
        dscr = 30e-7;       % Approx space charge region thickness
        pscr = 30;          % No of points in space charge region
        pepe = 20;          % electrode interface points
        te = 10e-7;         % electrode interface thickness
        dmin = 1e-7;        % start value for log meshes
        deltax = 1e-7;
        
        % Define spatial cordinate system
        % m=0 cartesian
        % m=1 cylindrical polar coordinates
        % m=2 spherical polar coordinates
        m = 0;
        
        %%%%%%% GENERAL CONTROL PARAMETERS %%%%%%%%%%
        OC = 0;                 % Closed circuit = 0, Open Circuit = 1
        Int = 0;                % Bias Light intensity (Suns Eq.)
        pulseon = 0;            % Switch pulse on TPC or TPV
        Vapp = 0;               % Applied bias
        BC = 3;                 % Boundary Conditions. Must be set to one for first solution
        figson = 1;             % Toggle figures on/off
        meshx_figon = 0;        % Toggles x-mesh figures on/off
        mesht_figon = 0;        % Toggles t-mesh figures on/off
        side = 1;               % illumination side 1 = EE, 2 = SE
        calcJ = 0;              % Calculates Currents- slows down solving calcJ = 1, calculates DD currents at every position
        mobset = 1;             % Switch on/off electron hole mobility- MUST BE SET TO ZERO FOR INITIAL SOLUTION
        mobseti = 1;
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
        stack = {'PEDOT', 'MAPICl', 'PCBM'}
        
        %% Energy levels    
        EA = [0, -0.1, -0.2];
        IP = [-1.0, -1.1, -1.2];
        
        % PCBM: Sigma Aldrich https://www.sigmaaldrich.com/technical-documents/articles/materials-science/organic-electronics/pcbm-n-type-semiconductors.html
        
        %% Equilibrium Fermi energies - defines doping density
        E0 = [-0.7, -0.6, -0.5];
        %E0 = [-5.25, -4.6, -3.95];
        
        % Workfunction energies
        PhiA = -0.7;
        PhiC = -0.5;
        
        % Effective Density Of States
        % DIFFERENT eDOS IN DIFFERENT LAYERS AS YET UNTESTED!
        N0 = [1e19, 1e19, 1e19];
        % PEDOT eDOS: https://aip.scitation.org/doi/10.1063/1.4824104
        % MAPI eDOS: F. Brivio, K. T. Butler, A. Walsh and M. van Schilfgaarde, Phys. Rev. B, 2014, 89, 155204.
        % PCBM eDOS:
                
        %%%%% MOBILE ION DEFECT DENSITY %%%%%
        Nion = [0, 1e18, 0];                     % [cm-3] A. Walsh, D. O. Scanlon, S. Chen, X. G. Gong and S.-H. Wei, Angewandte Chemie, 2015, 127, 1811.
        DOSion = [0, 1.21e22, 0];                % [cm-3] max density of iodide sites- P. Calado thesis
        
        % Mobilities
        mue = [1, 1, 1];         % electron mobility [cm2V-1s-1]
        muh = [1, 1, 1];         % hole mobility [cm2V-1s-1]
        muion = [0, 1e-10, 0];                   % ion mobility [cm2V-1s-1]
        % PTPD h+ mobility: https://pubs.rsc.org/en/content/articlehtml/2014/ra/c4ra05564k
        % PEDOT mue = 0.01 cm2V-1s-1 https://aip.scitation.org/doi/10.1063/1.4824104
        % TiO2 mue = 0.09?Bak2008
        % Spiro muh = 0.02 cm2V-1s-1 Hawash2018
        
        % Dielectric constants
        epp = [4,4,4];
        % TiO2?Wypych2014
        
        %%% Generation %%%%%
        G0 = [0, 2.6409e+21, 0];        % Uniform generation rate @ 1 Sun      
        
        %%%%%%% RECOMBINATION %%%%%%%%%%%
        % Radiative recombination, U = k(np - ni^2)
        krad = [1e-12, 1e-12, 1e-12];  % [cm3 s-1] Radiative Recombination coefficient
        % p-type/pi-interface/intrinsic/in-interface/n-type
                
        % SRH trap energies- currently set to mid gap - always set
        % respective to energy levels to avoid conflicts
        % U = (np-ni^2)/(taun(p+pt) +taup(n+nt))           
        Et =[-0.5, -0.5, -0.5];
    
        taun = [1e-9, 1e-6, 1e-12];           % [s] SRH time constant for electrons
        taup = [1e-9, 1e-6, 1e-12];           % [s] SRH time constant for holes
        
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
        
        dev;
        xx;
                
    end
    
    
    %%  Properties whose values depend on other properties (see 'get' methods).
    properties (Dependent)
        
        dcum
        dEAdx
        dIPdx
        dN0dx
        Eg
        Eif
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
            
            params.dev  = pc.builddev(params);
            params.xx = pc.xmeshini(params);
            
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
            
            % warn if trap energies are outside of band gpa energies
            for i = 1:length(params.Et)
                if params.Et(i) >= params.EA(2) || params.Et(i) <= params.IP(2)
                    msg = 'Trap energies must exist within layer 2 band gap.';
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
                               
        % Cumulative thickness
        function value = get.dcum(params)
                value = cumsum(params.d);                % cumulative thickness     
        end
        
        % Band gap energies
        function value = get.Eg(params)
            
            value = params.EA - params.IP;
            
        end
        
        % Vbi based on difference in boundary workfunctions
        function value = get.Vbi(params)
            
            value = params.PhiC - params.PhiA;
            
        end
                
        % Intrinsic Fermi Energies (Botzmann)
        function value = get.Eif(params)
            
            value = [0.5*((params.EA(1)+params.IP(1))+params.kB*params.T*log(params.N0(1)/params.N0(1))),...
                0.5*((params.EA(2)+params.IP(2))+params.kB*params.T*log(params.N0(2)/params.N0(2))),...
                0.5*((params.EA(3)+params.IP(3))+params.kB*params.T*log(params.N0(3)/params.N0(3)))];
            
        end
        
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
            
            value = [0, 0, F.fdn(params.N0(3), params.EA(3), params.E0(3), params.T)];

        end
        % Donor densities
        function value = get.NA(params)
            
            value = [F.fdp(params.N0(1), params.IP(1), params.E0(1), params.T), 0, 0];
            
        end
        
        % Intrinsic carrier densities - Boltzmann for now
        function value = get.ni(params)
            
            value = [params.N0(1)*(exp(-params.Eg(1)/(2*params.kB*params.T))),...
                params.N0(2)*(exp(-params.Eg(2)/(2*params.kB*params.T))),...
                params.N0(3)*(exp(-params.Eg(3)/(2*params.kB*params.T)))];
            
        end
        % Equilibrium electron densities
        function value = get.n0(params)
            
            value = [F.fdn(params.N0(1), params.EA(1), params.E0(1), params.T),...
                F.fdn(params.N0(2), params.EA(2), params.E0(2), params.T),...
                F.fdn(params.N0(3), params.EA(3), params.E0(3), params.T)];           
        end
                     
        function value = get.p0(params)
            
            value = [F.fdp(params.N0(1), params.IP(1), params.E0(1), params.T),...
                F.fdp(params.N0(2), params.IP(2), params.E0(2), params.T),...
                F.fdp(params.N0(3), params.IP(3), params.E0(3), params.T)];          % Background density holes in htl/p-type
            
        end
                
        % Boundary electron and hole densities - based on the workfunction
        % of the electrodes
        function value = get.nleft(params)
            
            value = F.fdn(params.N0(1), params.EA(1), params.PhiA, params.T);
        
        end
        
                % Boundary electron and hole densities
        function value = get.nright(params)
            
            value = F.fdn(params.N0(3), params.EA(3), params.PhiC, params.T);
        
        end
        
                % Boundary electron and hole densities
        function value = get.pleft(params)
            
            value = F.fdp(params.N0(1), params.IP(1), params.PhiA, params.T);
        
        end
        
        function value = get.pright(params)
            
            value = F.fdp(params.N0(3), params.IP(3), params.PhiC, params.T);
        
        end
        
        %% SRH parameters
        % nt_p - Density of CB electrons when Fermi level at trap state energy
        function value = get.nt(params)
            
            value = F.fdn(params.N0, params.EA, params.Et, params.T);
                        
        end
        
        % pt_p - Density of VB holes when Fermi level at trap state energy
        function value = get.pt(params)
            
            value =  F.fdp(params.N0, params.IP, params.Et, params.T);

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
    
    methods (Static)
        
        function xx = xmeshini(params)
            
            xx = meshgen_x(params);
            
        end
        
        % EA array
        function dev = builddev(params)
            % This function builds the properties for the device as
            % concatenated arrays such that each property can be called for
            % each point including grading at interfaces. For future
            % versions a choice of functions defining how the properties change
            % at the interfaces is intended. For the time being the
            % properties simply change linearly.
            xx = pc.xmeshini(params);
            
            dev.EA = zeros(1, length(xx));
            dev.IP = zeros(1, length(xx)); 
            dev.mue = zeros(1, length(xx));
            dev.muh = zeros(1, length(xx));
            dev.muion = zeros(1, length(xx));
            dev.NA = zeros(1, length(xx));
            dev.ND = zeros(1, length(xx));
            dev.N0 = zeros(1, length(xx));
            dev.Nion = zeros(1, length(xx));
            dev.ni = zeros(1, length(xx));
            dev.n0 = zeros(1, length(xx));
            dev.p0 = zeros(1, length(xx));
            dev.DOSion = zeros(1, length(xx));
            dev.epp = zeros(1, length(xx));
            dev.krad = zeros(1, length(xx));
            dev.gradEA = zeros(1, length(xx));
            dev.gradIP = zeros(1, length(xx));
            dev.gradN0 = zeros(1, length(xx));
            dev.E0 = zeros(1, length(xx));
            dev.G0 = zeros(1, length(xx));
            
            dcum0 = [0,params.dcum];
          for i =1:length(params.dcum)
                    if i == 1
                        aa = 0;
                    else
                        aa = 1;
                    end
                    
                    if i == length(params.dcum)
                        bb = 0;
                    else
                        bb = 1;
                    end
                    
                    for j = 1:length(xx)                     
                        if xx(j) >= dcum0(i) + aa*params.dint && xx(j) <= dcum0(i+1) - bb*params.dint                         
                            dev.EA(j) = params.EA(i); 
                            dev.IP(j) = params.IP(i);
                            dev.mue(j) = params.mue(i);
                            dev.muh(j) = params.muh(i);
                            dev.muion(j) = params.muion(i);
                            dev.N0(j) = params.N0(i);
                            dev.NA(j) = params.NA(i);
                            dev.ND(j) = params.ND(i);
                            dev.epp(j) = params.epp(i);
                            dev.ni(j) = params.ni(i);
                            dev.Nion(j) = params.Nion(i);
                            dev.DOSion(j) = params.DOSion(i);
                            dev.krad(j) = params.krad(i);
                            dev.n0(j) = params.n0(i);
                            dev.p0(j) = params.p0(i);
                            dev.E0(j) = params.E0(i);
                            dev.G0(j) = params.G0(i);
                            dev.gradEA(j) = 0;
                            dev.gradIP(j) = 0;
                            dev.gradN0(j) = 0;
                            
                        elseif xx(j) > dcum0(i+1) - bb*params.dint && xx(j) < dcum0(i+2) + bb*params.dint             
                            xprime = xx(j)-(dcum0(i+1) - bb*params.dint);
                            
                            dEAdxprime = (params.EA(i+1)-params.EA(i))/(2*params.dint);
                            dev.EA(j) = params.EA(i) + xprime*dEAdxprime;
                            dev.gradEA(j) = dEAdxprime;
                            
                            dIPdxprime = (params.IP(i+1)-params.IP(i))/(2*params.dint);
                            dev.IP(j) = params.IP(i) + xprime*dIPdxprime;
                            dev.gradIP(j) = dIPdxprime;
                                                        
                            % Mobilities
                            dmuedx = (params.mue(i+1)-params.mue(i))/(2*params.dint);
                            dev.mue(j) = params.mue(i) + xprime*dmuedx;
                            
                            dmuhdx = (params.muh(i+1)-params.muh(i))/(2*params.dint);
                            dev.muh(j) = params.muh(i) + xprime*dmuhdx;                           
                            
                            dmuiondx = (params.muion(i+1)-params.muion(i))/(2*params.dint);
                            dev.muion(j) = params.muion(i) + xprime*dmuiondx;
                            
                            dmuhdx = (params.muh(i+1)-params.muh(i))/(2*params.dint);
                            dev.muh(j) = params.muh(i) + xprime*dmuhdx;                                
                            
                            % effective density of states
                            dN0dx = (params.N0(i+1)-params.N0(i))/(2*params.dint);
                            dev.N0(j) = params.N0(i) + xprime*dN0dx;  
                            dev.gradN0(j) = dN0dx;
                            
                            % Doping densities
                            dNAdx = (params.NA(i+1)-params.NA(i))/(2*params.dint);
                            dev.NA(j) = params.NA(i) + xprime*dNAdx;   
                            
                            dNDdx = (params.ND(i+1)-params.ND(i))/(2*params.dint);
                            dev.ND(j) = params.ND(i) + xprime*dNDdx;
                            
                            % Dielectric constants
                            deppdx = (params.epp(i+1)-params.epp(i))/(2*params.dint);
                            dev.epp(j) = params.epp(i) + xprime*deppdx;
                            
                            % Intrinsic carrier densities
                            dnidx = (params.ni(i+1)-params.ni(i))/(2*params.dint);
                            dev.ni(j) = params.ni(i) + xprime*dnidx;
                            
                            % Equilibrium carrier densities
                            dn0dx = (params.n0(i+1)-params.n0(i))/(2*params.dint);
                            dev.n0(j) = params.n0(i) + xprime*dn0dx;
                                                                                    
                            % Equilibrium carrier densities
                            dp0dx = (params.p0(i+1)-params.p0(i))/(2*params.dint);
                            dev.p0(j) = params.p0(i) + xprime*dp0dx;
                            
                            % Equilibrium Fermi energy
                            dE0dx = (params.E0(i+1)-params.E0(i))/(2*params.dint);
                            dev.E0(j) = params.E0(i) + xprime*dE0dx;
                            
                            % Uniform generation rate
                            dG0dx = (params.G0(i+1)-params.G0(i))/(2*params.dint);
                            dev.G0(j) = params.G0(i) + xprime*dG0dx;
                            
                            % Static ion background density
                            dNiondx = (params.Nion(i+1)-params.Nion(i))/(2*params.dint);
                            dev.Nion(j) = params.Nion(i) + xprime*dNiondx;
                            
                            % Ion density of states
                            dDOSiondx = (params.DOSion(i+1)-params.DOSion(i))/(2*params.dint);
                            dev.DOSion(j) = params.DOSion(i) + xprime*dDOSiondx;
                            
                            %recombination
                            dkraddx = (params.krad(i+1)-params.krad(i))/(2*params.dint);
                            dev.krad(j) = params.krad(i) + xprime*dkraddx;
                           
                        end  
                    end                                   
                    
                end
                % calculate the gradients
                %dev.gradEA = gradient(dev.EA, xx);
                %dev.gradIP = gradient(dev.IP, xx);
                %dev.gradN0 = gradient(dev.N0, xx);
        end
        
    end
end
