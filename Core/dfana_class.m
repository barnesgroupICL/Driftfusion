classdef dfana_class
    
    methods (Static)
        
        function [u,t,x,par,dev,n,p,a,V] = splitsol(sol)
            % splits solution into useful outputs
            u = sol.u;
            t = sol.t;
            x = sol.x;
            par = sol.par;
            dev = par.dev;
            
            % split the solution into its component parts (e.g. electrons, holes and efield)
            n = u(:,:,1);
            p = u(:,:,2);
            a = u(:,:,3);
            V = u(:,:,4);
        end
        
        function [Ecb, Evb, Efn, Efp] = QFLs(sol)
            % u is the solution structure
            % Simple structure names
            [u,t,x,par,dev,n,p,a,V] = dfana_class.splitsol(sol);
                       
            % Create 2D matrices for multiplication with solutions
            EAmat = repmat(par.dev.EA, length(t), 1);
            IPmat = repmat(par.dev.IP, length(t), 1);
            muemat = repmat(par.dev.mue, length(t), 1);
            muhmat = repmat(par.dev.muh, length(t), 1);
            muionmat = repmat(par.dev.muion, length(t), 1);
            NAmat = repmat(par.dev.NA, length(t), 1);
            NDmat = repmat(par.dev.ND, length(t), 1);
            Ncmat = repmat(par.dev.Nc, length(t), 1);
            Nvmat = repmat(par.dev.Nv, length(t), 1);
            Nionmat = repmat(par.dev.Nion, length(t), 1);
            eppmat = repmat(par.dev.epp, length(t), 1);
            nimat = repmat(par.dev.ni, length(t), 1);
            kradmat = repmat(par.dev.krad, length(t), 1);
            taunmat = repmat(par.dev.taun, length(t), 1);
            taupmat = repmat(par.dev.taup, length(t), 1);
            ntmat = repmat(par.dev.nt, length(t), 1);
            ptmat = repmat(par.dev.pt, length(t), 1);
            
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
            [u,t,x,par,dev,n,p,a,V] = dfana_class.splitsol(sol);
            
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
            end
            
            dndtInt = trapz(x, dndt, 2);
            dpdtInt = trapz(x, dpdt, 2);
            dadtInt = trapz(x, dadt, 2);
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
            
            djndx = -(dndt - g + U);    % Not certain about the sign here
            djpdx = -(dpdt - g + U);
            djadx = -dadt;
            
            % Integrate across the device to get delta fluxes at all positions
            deltajn = cumtrapz(x, djndx, 2);
            deltajp = cumtrapz(x, djpdx, 2);
            deltaja = cumtrapz(x, djadx, 2);
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
            j.p = jp_l + deltajp;
            j.a = 0 + deltaja;
            
            % displacement flux
            j.disp = zeros(length(t), length(x));
            [FV, Frho] = dfana_class.calcF(sol);
            
            for i = 1:length(x)
                j.disp(:,i) = par.epp0.*eppmat(:,i).*(gradient(Frho(:,i), t));
            end
            
            J.n = -j.n*par.e;
            J.p = j.p*par.e;
            J.a = j.a*par.e;
            J.disp = j.disp*abs(par.e);
            
            % Total current
            J.tot = J.n + J.p + J.a + J.disp;
        end
        
        function [FV, Frho] = calcF(sol)
            % Electric field caculation
            % FV = Field calculated from the gradient of the potential
            % Frho = Field calculated from integrated space charge density  
            [u,t,x,par,dev,n,p,a,V] = dfana_class.splitsol(sol);
            rho = dfana_class.calcrho(sol);
            eppmat = repmat(dev.epp, length(t), 1);
            
            for i=1:length(t)
                FV(i,:) = -gradient(V(i, :), x);                      % Electric field calculated from V
            end
            
            Frho = cumtrapz(x, rho./(eppmat.*par.epp0), 2) + FV(:,1);
            
        end
        
        
        function rho = calcrho(sol)
            % Calculates the space charge density
            [u,t,x,par,dev,n,p,a,V] = dfana_class.splitsol(sol);
            
            NAmat = repmat(dev.NA, length(t), 1);
            NDmat = repmat(dev.ND, length(t), 1);
            Nionmat = repmat(dev.Nion, length(t), 1);
            
            % charge density
            rho = -n + p + a -NAmat + NDmat - Nionmat;
            
        end
        
        
        
        
    end
    
end
