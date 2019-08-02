classdef dfana
    % Analysis class
    methods (Static)
        
        function [u,t,x,par,dev,n,p,a,c,V] = splitsol(sol)
            % splits solution into useful outputs
            u = sol.u;
            t = sol.t;
            x = sol.x;
            par = sol.par;
            dev = par.dev;
            
            % split the solution into its component parts (e.g. electrons, holes and efield)
            n = u(:,:,1);
            p = u(:,:,2);
            c = u(:,:,3);
            V = u(:,:,4);
            if par.N_ionic_species == 2
                a = u(:,:,5);
            else
                a = repmat(dev.Ncat, length(t), 1);
            end
        end
        
        function [Ecb, Evb, Efn, Efp] = QFLs(sol)
            % u is the solution structure
            % Simple structure names
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            % Create 2D matrices for multiplication with solutions
            EAmat = repmat(dev.EA, length(t), 1);
            IPmat = repmat(dev.IP, length(t), 1);
            NAmat = repmat(dev.NA, length(t), 1);
            NDmat = repmat(dev.ND, length(t), 1);
            Ncmat = repmat(dev.Nc, length(t), 1);
            Nvmat = repmat(dev.Nv, length(t), 1);
            Nionmat = repmat(dev.Nion, length(t), 1);
            Ncatmat = repmat(dev.Ncat, length(t), 1);
            eppmat = repmat(dev.epp, length(t), 1);
            nimat = repmat(dev.ni, length(t), 1);
            
            Ecb = EAmat-V;                                 % Conduction band potential
            Evb = IPmat-V;                                 % Valence band potential
            
            Efn = zeros(size(n,1), size(n,2));
            Efp = zeros(size(n,1), size(n,2));
            
            if par.stats == 'Fermi'
                
                for i = 1:size(n,1)           % time
                    for j = 1:size(n,2)       % position
                        Efn(i,j) = F.Efn_fd_fun(n(i,j), dev.Efn(j,:),  dev.n_fd(j,:));
                        Efp(i,j) = F.Efp_fd_fun(p(i,j), dev.Efp(j,:),  dev.p_fd(j,:));
                    end
                end
                Efn = Efn-V;
                Efp = Efp-V;
                
            elseif par.stats == 'Boltz'
                Efn = real(Ecb+(par.kB*par.T/par.q)*log(n./Ncmat));        % Electron quasi-Fermi level
                Efp = real(Evb-(par.kB*par.T/par.q)*log(p./Nvmat));        % Hole quasi-Fermi level
            end
            
        end
        
        function [j, J] = calcJ(sol)
            % Current, J and flux, j calculation from continuity equations
            
            % obtain SOL components for easy referencing
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            % Read in generation profil
            if par.OM == 1
                gx = sol.gx;
            end
            
            % Property matrices
            eppmat = repmat(dev.epp, length(t), 1);
            nimat = repmat(dev.ni, length(t), 1);
            kradmat = repmat(dev.krad, length(t), 1);
            taunmat = repmat(dev.taun, length(t), 1);
            taupmat = repmat(dev.taup, length(t), 1);
            ntmat = repmat(dev.nt, length(t), 1);
            ptmat = repmat(dev.pt, length(t), 1);
            
            for i = 1:size(n, 2)
                dndt(:,i) = gradient(n(:,i), t);
                dpdt(:,i) = gradient(p(:,i), t);
                dadt(:,i) = gradient(a(:,i), t);
                dcdt(:,i) = gradient(c(:,i), t);
            end
            
            dndtInt = trapz(x, dndt, 2);
            dpdtInt = trapz(x, dpdt, 2);
            dadtInt = trapz(x, dadt, 2);
            dcdtInt = trapz(x, dcdt, 2);
            % Recombination
            Ubtb = kradmat.*(n.*p - nimat.^2);
            
            Usrh = ((n.*p - nimat.^2)./((taunmat.*(p+ptmat)) + (taupmat.*(n+ntmat))));
            
            U = Ubtb + Usrh;
            
            Usrhnogen = ((n.*p)./((taunmat.*(p+ptmat)) + (taupmat.*(n+ntmat))));
            
            % Uniform Generation
            switch par.OM
                
                % Uniform generation
                case 0
                    
                    g = par.Int*dev.G0;
                    
                case 1
                    
                    gxAM15 = par.Int*repmat(gx.AM15', length(t), 1);
                    
                    if par.pulseon == 1
                        
                        las = repmat(gx.las', length(t), 1);
                        pulset = ones(length(x), length(t))*diag(t >= par.pulsestart & t < par.pulselen + par.pulsestart);
                        pulset = pulset';
                        gxlas = las.*pulset;
                        
                    else
                        gxlas = 0;
                        
                    end
                    
                    g = gxAM15 + gxlas;
                    
                case 2
                    % Transfer Matrix
                    if par.Int == 0
                        
                        g = 0;
                        
                    else
                        
                        g = par.Int*interp1(par.genspace, solstruct.Gx1S, (x-par.dcum(1)));
                        
                    end
                    
            end
            
            djndx = (dndt + g - U);    % Not certain about the sign here
            djpdx = (dpdt + g - U);
            djadx = dadt;
            djcdx = dcdt;
            
            % Integrate across the device to get delta fluxes at all positions
            deltajn = cumtrapz(x, djndx, 2);
            deltajp = cumtrapz(x, djpdx, 2);
            deltaja = cumtrapz(x, djadx, 2);
            deltajc = cumtrapz(x, djcdx, 2);
            %% Currents from the boundaries
            switch par.BC
                case 0
                    jn_l = 0;
                    jp_l = 0;
                    jn_r = 0;
                    jp_r = 0;
                    % Blocking contacts
                case 1
                    % Setting jp_l = djpdx(end) ensures that jp_r = 0;
                    jn_l = 0;
                    jp_l = -deltajp(:, end);
                    
                    jn_r = deltajn(:, end);
                    jp_r = 0;
                    
                case 2
                    
                    jn_l = -par.sn_l*(n(:, 1) - par.nleft);
                    jp_l = -deltajp(:, end) + par.sp_r*(p(:, end) - par.pright);
                    
                    jn_r = deltajn(:, end) - par.sn_l*(n(:, 1) - par.nleft);
                    jp_r = par.sp_r*(p(:, end) - par.pright);
                    
                case 3
                    
                    jn_l = -par.sn_l*(n(:, 1) - par.nleft);
                    jp_l = -par.sp_l*(p(:, 1) - par.pleft);
                    
                    jn_r = par.sn_r*(n(:, end) - par.nright);
                    jp_r = par.sp_r*(p(:, end) - par.pright);
                    
            end
            
            % Calculate total electron and hole currents from fluxes
            j.n = jn_l + deltajn;
            j.p = jp_r + (deltajp-deltajp(end));
      
            j.a = 0 + deltaja;
            j.c = 0 + deltajc;
            % displacement flux
            j.disp = zeros(length(t), length(x));
            [FV, Frho] = dfana.calcF(sol);
            
            for i = 1:length(x)
                j.disp(:,i) = par.epp0.*eppmat(:,i).*(gradient(Frho(:,i), t));
            end
            
            J.n = j.n*-par.e;
            J.p = j.p*par.e;
            J.a = j.a*-par.e;
            J.c = j.c*par.e;
            J.disp = j.disp*abs(par.e);
            
            % Total current
            J.tot = J.n + J.p + J.a + J.c + J.disp;
        end
        
        function U = calcU(sol)
            % obtain SOL components for easy referencing
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            % Read in generation profil
            if par.OM == 1
                gx = sol.gx;
            end
            
            % Property matrices
            eppmat = repmat(dev.epp, length(t), 1);
            nimat = repmat(dev.ni, length(t), 1);
            kradmat = repmat(dev.krad, length(t), 1);
            taunmat = repmat(dev.taun, length(t), 1);
            taupmat = repmat(dev.taup, length(t), 1);
            ntmat = repmat(dev.nt, length(t), 1);
            ptmat = repmat(dev.pt, length(t), 1);
            
            for i = 1:size(n, 2)
                dndt(:,i) = gradient(n(:,i), t);
                dpdt(:,i) = gradient(p(:,i), t);
                dadt(:,i) = gradient(a(:,i), t);
                dcdt(:,i) = gradient(c(:,i), t);
            end
            
            dndtInt = trapz(x, dndt, 2);
            dpdtInt = trapz(x, dpdt, 2);
            dadtInt = trapz(x, dadt, 2);
            dcdtInt = trapz(x, dcdt, 2);
            % Recombination
            U.btb = kradmat.*(n.*p - nimat.^2);
            
            U.srh = ((n.*p - nimat.^2)./((taunmat.*(p+ptmat)) + (taupmat.*(n+ntmat))));
            
            U.tot = U.btb + U.srh;
        end
        
        function Jdd = Jddxt(sol)
            % obtain SOL components for easy referencing
            [u,t,xmesh,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            xhalfmesh = getxihalf(sol);
            devihalf = getdevihalf(par);
            
            % Property matrices
            eppmat = devihalf.epp;
            mue_mat = devihalf.mue;
            muh_mat = devihalf.muh;
            mue_cat = devihalf.mucat;
            muh_ion = devihalf.muion;
            gradEA_mat = devihalf.gradEA;
            gradIP_mat = devihalf.gradIP;
            gradNc_mat = devihalf.gradNc;
            gradNv_mat = devihalf.gradNv;
            Nc_mat = devihalf.Nc;
            Nv_mat = devihalf.Nv;
            % Calculates drift and diffusion currents at every point and all times -
            % NOTE: UNRELIABLE FOR TOTAL CURRENT as errors in the calculation of the
            % spatial gradients mean that the currents do not cancel properly
            for i = 1:length(t)
                   
                [nloc(i,:),dndx(i,:)] = pdeval(0,xmesh,n(i,:),xhalfmesh);
                [ploc(i,:),dpdx(i,:)] = pdeval(0,xmesh,p(i,:),xhalfmesh);
                [aloc(i,:),dadx(i,:)] = pdeval(0,xmesh,a(i,:),xhalfmesh);
                [cloc(i,:),dcdx(i,:)] = pdeval(0,xmesh,c(i,:),xhalfmesh);
                [Vloc(i,:),dVdx(i,:)] = pdeval(0,xmesh,V(i,:),xhalfmesh);
    
                % Diffusion coefficients
                if par.stats == 'Fermi'
                    for jj = 1:length(x)
                        Dn(i,jj) = F.D(nloc(i,jj), devihalf.Dnfun(jj,:), devihalf.n_fd(jj,:));
                        Dp(i,jj) = F.D(ploc(i,jj), devihalf.Dpfun(jj,:), devihalf.p_fd(jj,:));
                    end
                end
            end
            
            if par.stats == 'Boltz'
                Dn_mat = mue_mat*par.kB*par.T;
                Dp_mat = muh_mat*par.kB*par.T;
            end
            
            % Particle currents
%             Jdd.ndiff = -Dn_mat.*(dndx-((nloc./Nc_mat).*gradNc_mat)).*-par.e;
%             Jdd.ndrift = mue_mat.*nloc.*(dVdx-gradEA_mat)*-par.e;
%             
%             Jdd.pdiff = -Dp_mat.*(dpdx-((ploc./Nv_mat).*gradNv_mat)).*par.e;
%             Jdd.pdrift = muh_mat.*ploc.*(-dVdx+gradIP_mat)*par.e;

            Jdd.ndiff = -Dn_mat.*(dndx-((nloc./Nc_mat).*gradNc_mat)).*-par.e;
            Jdd.ndrift = mue_mat.*nloc.*(dVdx-gradEA_mat)*-par.e;
            
            Jdd.pdiff = -Dp_mat.*(dpdx-((ploc./Nv_mat).*gradNv_mat)).*par.e;
            Jdd.pdrift = muh_mat.*ploc.*(-dVdx+gradIP_mat)*par.e;

            Jdd.cdiff = -devihalf.mucat.*par.kB*par.T.*dcdx*par.e;
            Jdd.cdrift = devihalf.mucat.*cloc.*-dVdx*par.e;
            
            if par.N_ionic_species == 1
                Jdd.adiff = zeros(length(t), length(xhalfmesh));
                Jdd.adrift = zeros(length(t), length(xhalfmesh));
            elseif par.No_ionic_species == 2
                Jdd.adiff = -devihalf.muion.*par.kB*par.T.*dadx.*-par.e;
                Jdd.adrift = devihalf.muion.*aloc.*dVdx.*-par.e;
            end
            
            Jdd.n = Jdd.ndrift + Jdd.ndiff;
            Jdd.p = Jdd.pdrift + Jdd.pdiff;
            Jdd.a = Jdd.adrift + Jdd.adiff;
            Jdd.c = Jdd.cdrift + Jdd.cdiff;
            
            % Displacement current
            %[FV, Frho] = dfana.calcF(sol);
            for i = 1:length(xhalfmesh)
                j.disp(:,i) = par.epp0.*eppmat(:,i).*(gradient(dVdx(:,i), t));
            end
            
            Jdd.disp = j.disp*abs(par.e);
            
            Jdd.tot = Jdd.n + Jdd.p + Jdd.a + Jdd.c + Jdd.disp;
        end
        
        function [FV, Frho] = calcF(sol)
            % Electric field caculation
            % FV = Field calculated from the gradient of the potential
            % Frho = Field calculated from integrated space charge density
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            rho = dfana.calcrho(sol);
            eppmat = repmat(dev.epp, length(t), 1);
            
            for i=1:length(t)
                FV(i,:) = -gradient(V(i, :), x);                      % Electric field calculated from V
            end
            
            Frho = cumtrapz(x, rho./(eppmat.*par.epp0), 2) + FV(:,1);
            
        end
        
        function rho = calcrho(sol)
            % Calculates the space charge density
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            NAmat = repmat(dev.NA, length(t), 1);
            NDmat = repmat(dev.ND, length(t), 1);
            
            % charge density
            rho = -n + p - a + c - NAmat + NDmat;
        end
        
        
        function Vapp = calcVapp(sol)
            par = sol.par;
            if par.JV == 2
                Vapp = par.Vapp_func(par.Vapp_params, sol.t);
            else
                Vapp = -(sol.u(:,end,4)-sol.u(:,1,4)-sol.par.Vbi);
            end
        end
        
        function stats = JVstats(JVsol)
            % A function to pull statistics from a JV sweep using DOJV
            % JVsol - a solution from DOJV
            if isfield(JVsol, 'ill')
                if isfield(JVsol.ill, 'f')
                    Vapp = dfana.calcVapp(JVsol.ill.f);
                    [j,J] = dfana.calcJ(JVsol.ill.f);
                    try
                        stats.Jsc_f = interp1(Vapp, J.tot(:, end), 0);
                    catch
                        warning('No Jsc available- Vapp must pass through 0 to obtain Jsc')
                        stats.Jsc_f = 0;
                    end
                    
                    try
                        stats.Voc_f = interp1(J.tot(:, end), Vapp, 0);
                    catch
                        warning('No Voc available- try increasing applied voltage range')
                        stats.Voc_f = 0;
                    end
                    
                    if stats.Jsc_f ~= 0 && stats.Voc_f ~= 0
                        pow_f = J.tot(:,end).*Vapp;
                        stats.mpp_f = min(pow_f);
                        stats.mppV_f = Vapp(pow_f == stats.mpp_f);
                        stats.FF_f = stats.mpp_f/(stats.Jsc_f*stats.Voc_f);
                    end
                    
                    %% Hysteresis Index
                    A_f = abs(trapz(Vapp(Vapp >=0 & Vapp <= stats.Voc_f), J.tot(Vapp >= 0 & Vapp <= stats.Voc_f, end)));
                    
                else
                    stats.Jsc_f = nan;
                    stats.Voc_f = nan;
                    stats.mpp_f = nan;
                    stats.FF_f = nan;
                end
                
                if isfield(JVsol.ill, 'r')
                    Vapp = dfana.calcVapp(JVsol.ill.r);
                    [j,J] = dfana.calcJ(JVsol.ill.r);
                    try
                        stats.Jsc_r = interp1(Vapp, J.tot(:, end), 0);
                    catch
                        warning('No Jsc available- Vapp must pass through 0 to obtain Jsc')
                        stats.Jsc_r = 0;
                    end
                    
                    try
                        stats.Voc_r = interp1(J.tot(:, end), Vapp, 0);
                    catch
                        warning('No Voc available- try increasing applied voltage range')
                        stats.Voc_r = 0;
                    end
                    
                    if stats.Jsc_r ~= 0 && stats.Voc_r ~= 0
                        pow_r = J.tot(:,end).*Vapp;
                        stats.mpp_r = min(pow_r);
                        stats.mppV_r = Vapp(pow_r == stats.mpp_r);
                        stats.FF_r = stats.mpp_r/(stats.Jsc_r*stats.Voc_r);
                    end
                    
                    %% Hysteresis Factor
                    A_r = abs(trapz(Vapp(Vapp >=0 & Vapp <= stats.Voc_r), J.tot(Vapp >= 0 & Vapp <= stats.Voc_r, end)));
                    
                    %% Sign to identify inverted hysteresis
                    if A_r >= A_f
                        B = 1;
                    elseif A_r < A_f
                        B = -1;
                    end
                    
                    stats.HF = B*abs((A_r - A_f)/A_r);
                    
                else
                    stats.Jsc_r = NaN;
                    stats.Voc_r = NaN;
                    stats.mpp_r = NaN;
                    stats.FF_r = NaN;
                    stats.HF = NaN;
                end
            else
            end
        end
        
        function F = Ft(sol, ppos)
            % Field as a function of time at point position PPOS
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            for i = 1:length(t)
                F(i,:) = -gradient(V(i,:),x);
            end
            F = F(:,ppos);
        end
        
        function value = PLt(sol)
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            value = trapz(x,(n.*p),2);
        end
        
        function value = Voct(sol)
            %Get QFLs
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            value = Efn(:, end) - Efp(:, 1);
        end
        
        function deltaV = deltaVt(sol, p1, p2)
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            % Calculates the electorstatic potential difference as a function of time
            % between two points P1 and P2
            deltaV = V(:,p1) - V(:,p2);
        end
        
        function sigma = calcsigma(sol)
            % calculates the integrated space charge density
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            rho = dfana.calcrho(sol);
            sigma = trapz(x, rho, 2);
        end
    end
end
