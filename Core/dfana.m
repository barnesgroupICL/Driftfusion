classdef dfana
    % DRIFTFUSION analysis class- contains multiple methods for calculating
    % outputs using the solution obtained from DF.
    %
    %% LICENSE
    % Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
    % Imperial College London
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU Affero General Public License as published
    % by the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    %% Start code
    methods (Static)
        function [u,t,x,par,dev,n,p,a,c,V] = splitsol(sol)
            % splits solution into useful outputs
            u = sol.u;
            t = sol.t(1:size(u,1));
            x = sol.x;
            par = sol.par;
            dev = par.dev;

            % split the solution into its component parts (e.g. electrons, holes and efield)
            V = u(:,:,1);
            n = u(:,:,2);
            p = u(:,:,3);

            if par.N_ionic_species == 1
                c = u(:,:,4);
                a = repmat(dev.Nani, length(t), 1);
            elseif par.N_ionic_species == 2
                c = u(:,:,4);
                a = u(:,:,5);
            else
                c = repmat(dev.Ncat, length(t), 1);
                a = repmat(dev.Nani, length(t), 1);
            end
        end

        function [Ecb, Evb, Efn, Efp] = QFLs(sol)
            % u is the solution structure
            % Simple structure names
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            Ecb = dev.EA-V;                                 % Conduction band potential
            Evb = dev.IP-V;                                 % Valence band potential
            Efn = zeros(size(n,1), size(n,2));
            Efp = zeros(size(n,1), size(n,2));

            if par.prob_distro_function == 'Fermi'

                for i = 1:size(n,1)           % time
                    for j = 1:size(n,2)       % position
                        Efn(i,j) = distro_fun.Efn_fd_fun(n(i,j), dev.Efn(j,:),  dev.n_fd(j,:));
                        Efp(i,j) = distro_fun.Efp_fd_fun(p(i,j), dev.Efp(j,:),  dev.p_fd(j,:));
                    end
                end
                Efn = Efn-V;
                Efp = Efp-V;

            elseif par.prob_distro_function == 'Boltz'
                Efn = real(Ecb+(par.kB*par.T/par.q)*log(n./dev.Nc));        % Electron quasi-Fermi level
                Efp = real(Evb-(par.kB*par.T/par.q)*log(p./dev.Nv));        % Hole quasi-Fermi level
            end

        end

        function [Ecb, Evb, Efn, Efp] = QFL_ihalf(sol)
            % Calculates the QFLs on the i_half xmesh
            % u is the solution structure
            % Simple structure names
            [u,t,xmesh,par,dev,n0,p0,a0,c0,V0] = dfana.splitsol(sol);


            n = getvarihalf(n0);
            p = getvarihalf(p0);
            V = getvarihalf(V0);
            dev_ihalf = par.dev_ihalf;
            Ecb_ihalf = dev_ihalf.EA-V;                                 % Conduction band potential
            Evb_ihalf = dev_ihalf.IP-V;                                 % Valence band potential

            if par.prob_distro_function == 'Fermi'

                for i = 1:size(n,1)           % time
                    for j = 1:size(n,2)       % position
                        Efn_ihalf(i,j) = distro_fun.Efn_fd_fun(n(i,j), dev.Efn(j,:),  dev.n_fd(j,:));
                        Efp_ihalf(i,j) = distro_fun.Efp_fd_fun(p(i,j), dev.Efp(j,:),  dev.p_fd(j,:));
                    end
                end
                Efn_ihalf = Efn_ihalf-V;
                Efp_ihalf = Efp_ihalf-V;

            elseif par.prob_distro_function == 'Boltz'
                Efn_ihalf = real(Ecb_ihalf+(par.kB*par.T/par.q)*log(n./dev_ihalf.Nc));        % Electron quasi-Fermi level
                Efp_ihalf = real(Evb_ihalf-(par.kB*par.T/par.q)*log(p./dev_ihalf.Nv));        % Hole quasi-Fermi level
            end
            Efn = Efn_ihalf;
            Efp = Efp_ihalf;
            Ecb = Ecb_ihalf;
            Evb = Evb_ihalf;
        end

        function [J, j, x] = calcJ(sol)
            % Current, J and flux, j calculation from continuity equations
            % obtain SOL components for easy referencing
            [u,t,xmesh,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            n_ihalf = getvarihalf(n);
            p_ihalf = getvarihalf(p);
            a_ihalf = getvarihalf(a);
            c_ihalf = getvarihalf(c);

            x = par.x_ihalf;
            [~,~,g] = dfana.calcg(sol);

            [~, dndt] = gradient(n_ihalf, x, t);
            [~, dpdt] = gradient(p_ihalf, x, t);
            [~, dadt] = gradient(a_ihalf, x, t);
            [~, dcdt] = gradient(c_ihalf, x, t);

            % Recombination
            r = dfana.calcr_ihalf(sol);

            djndx = -dndt + g - r.tot;
            djpdx = -dpdt + g - r.tot;
            djadx = -dadt;
            djcdx = -dcdt;

            deltajn = cumtrapz(x, djndx, 2);
            deltajp = cumtrapz(x, djpdx, 2);
            deltaja = cumtrapz(x, djadx, 2);
            deltajc = cumtrapz(x, djcdx, 2);

            %% Currents from the boundaries
            switch par.BC
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
            % Use the minority carrier flux as the boundary condition
            if par.pleft >= par.nleft && par.nright >= par.pright
                % p-type left boundary, n-type right boundary
                j.n = jn_l + deltajn;
                j.p = jp_r + (deltajp - deltajp(:,end));
            elseif par.nleft >= par.nright && par.pright >= par.nright
                % n-type left boundary, p-type right boundary
                j.n = jn_r + (deltajn-deltajn(:,end));
                j.p = jp_l + deltajp;
            elseif par.pleft >= par.nleft && par.pright >= par.nright...
                    || par.nleft >= par.pleft && par.nright >= par.pright
                % p-type both boundaries or n-type both boundaries
                j.n = jn_l + deltajn;
                j.p = jp_l + deltajp;
            end

            j.a = 0 + deltaja;
            j.c = 0 + deltajc;
            % displacement flux
            FV_ihalf = dfana.calcF_ihalf(sol);

            [~, FV_ihalf_dt] = gradient(FV_ihalf, x, t);
            j.disp = par.epp0.*par.dev_ihalf.epp.*FV_ihalf_dt;

            J.n = j.n*-par.e;
            J.p = j.p*par.e;
            J.a = j.a*-par.e;
            J.c = j.c*par.e;
            J.disp = j.disp*par.e;

            % Total current
            J.tot = J.n + J.p + J.a + J.c + J.disp;
        end

        function [g1, g2, g] = calcg(sol)
            [~,tmesh,~,par,~,~,~,~,~,~] = dfana.splitsol(sol);
            %% Generation function
            switch par.g1_fun_type
                case 'constant'
                    g1 = repmat(par.int1.*par.gx1, length(tmesh), 1);
                otherwise
                    g1_fun = fun_gen(par.g1_fun_type);
                    g1 = g1_fun(par.g1_fun_arg, tmesh')*par.gx1;
            end

            switch par.g2_fun_type
                case 'constant'
                    g2 = repmat(par.int2.*par.gx2, length(tmesh), 1);
                otherwise
                    g2_fun = fun_gen(par.g2_fun_type);
                    g2 = g2_fun(par.g2_fun_arg, tmesh')*par.gx2;
            end
            g = g1 + g2;
        end

        function r = calcr(sol)
            % Calculate the recombination rate on whole mesh
            % obtain SOL components for easy referencing
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            int_switch = repmat(dev.int_switch, length(t), 1);
            rec_zone = repmat(dev.int_switch, length(t), 1);
            bulk_switch = abs(int_switch-1);
            for i=1:length(t)
                dVdx(i,:) = pdeval(0, x, V(i,:), x);
            end
            xprime = dev.xprime;
            xprime_n = dev.xprime_n;
            xprime_p = dev.xprime_p;
            dint = dev.dint;
            alpha0 = dev.alpha0;
            beta0 = dev.beta0;
            alpha = par.q*dVdx./(par.kB*par.T) + alpha0;
            beta = par.q*dVdx./(par.kB*par.T) + beta0;

            % Band-to-band
            r.btb = dev.B.*(n.*p - dev.ni.^2);
            % Bulk SRH
            r.srh = bulk_switch.*((n.*p - dev.ni.^2)./((dev.taun.*(p+dev.pt)) + (dev.taup.*(n+dev.nt))));
            % Volumetric surface SRH
            if alpha < 0
                ns = n.*exp(-alpha.*xprime);
            elseif alpha >=0
                ns = n.*exp(-alpha.*(xprime-dint));
            end
            
            if beta < 0
                ps = p.*exp(-beta.*xprime);
            elseif beta >=0
                ps = p.*exp(-beta.*(xprime-dint));
            end
            
            r.vsr = rec_zone.*((ns.*ps - dev.ni.^2)./...
                ((dev.taun_vsr.*(ps + dev.pt)) + (dev.taup_vsr.*(ns + dev.nt))));
            % Total
            r.tot = r.btb + r.srh + r.vsr;
        end

        function [r, ns, ps] = calcr_ihalf(sol)
            % Calculate the recombination rate on i-half mesh
            % obtain SOL components for easy referencing
            [u,t,x,par,~,n,p,a,c,V] = dfana.splitsol(sol);
            dev = par.dev_ihalf;
            x_ihalf = par.x_ihalf;
            n_ihalf = getvarihalf(n);
            p_ihalf = getvarihalf(p);
            dVdx_ihalf = zeros(length(t), length(x_ihalf));
            for i=1:length(t)
                [~, dVdx_ihalf(i,:)] = pdeval(0, x, V(i,:), x_ihalf);
            end
            int_switch = repmat(dev.int_switch, length(t), 1);
            rec_zone = repmat(dev.rec_zone, length(t), 1);
            bulk_switch = abs(int_switch-1);
            xprime_n = dev.xprime_n;
            xprime_p = dev.xprime_p;
            xprime = dev.xprime;
            dint = dev.dint;
            
            alpha0 = dev.alpha0;
            beta0 = dev.beta0;

            alpha = abs(par.q*dVdx_ihalf./(par.kB*par.T) + alpha0);
            beta = abs(par.q*-dVdx_ihalf./(par.kB*par.T) + beta0);

            % Band-to-band
            r.btb = dev.B.*(n_ihalf.*p_ihalf - dev.ni.^2);
            % Bulk SRH
            r.srh = bulk_switch.*(n_ihalf.*p_ihalf - dev.ni.^2)...
                ./(dev.taun.*(p_ihalf+dev.pt) + dev.taup.*(n_ihalf+dev.nt));
            % Volumetric surface SRH
%             ns = n_ihalf.*exp(alpha.*xprime_n); % Projected electron surface density
%             ps = p_ihalf.*exp(beta.*xprime_p);  % Projected hole surface density
            if alpha < 0
                ns = n_ihalf.*exp(-alpha.*xprime);
            elseif alpha >=0
                ns = n_ihalf.*exp(-alpha.*(xprime-dint));
            end
            
            if beta < 0
                ps = p_ihalf.*exp(beta.*xprime);
            elseif beta >=0
                ps = p_ihalf.*exp(beta.*(xprime-dint));
            end
            r.vsr = rec_zone.*(ns.*ps - dev.ni.^2)...
                ./(dev.taun_vsr.*(ps + dev.pt) + dev.taup_vsr.*(ns + dev.nt));
            % Total
            r.tot = r.btb + r.srh + r.vsr;
        end

        function [jdd, Jdd, xout] = Jddxt(sol)
            % Calculates drift and diffusion currents at every point and all times -
            % NOTE: UNRELIABLE FOR TOTAL CURRENT as errors in the calculation of the
            % spatial gradients mean that the currents do not cancel properly
            % obtain SOL components for easy referencing
            [~,t,x,par,~,n,p,a,c,V] = dfana.splitsol(sol);
            xout = par.x_ihalf;
            dev = par.dev_ihalf;

            % Property matrices
            eppmat = dev.epp;
            mue_mat = dev.mue;
            muh_mat = dev.muh;
            mu_cat = dev.mucat;
            mu_ani = dev.muani;
            gradEA_mat = dev.gradEA;
            gradIP_mat = dev.gradIP;
            gradNc_mat = dev.gradNc;
            gradNv_mat = dev.gradNv;
            Nc_mat = dev.Nc;
            Nv_mat = dev.Nv;
            
            V_sub = zeros(length(t), length(xout));
            n_sub = zeros(length(t), length(xout));
            p_sub = zeros(length(t), length(xout));
            a_sub = zeros(length(t), length(xout));
            c_sub = zeros(length(t), length(xout));
            
            dVdx = zeros(length(t), length(xout));
            dndx = zeros(length(t), length(xout));
            dpdx = zeros(length(t), length(xout));
            dadx = zeros(length(t), length(xout));
            dcdx = zeros(length(t), length(xout));             
            
            %% Avoid PDEVAL for faster calculation
            % Obtain variables and gradients on sub-interval mesh
            for i = 1:length(t)
                V_sub(i,:) = 0.5*(V(i, 2:end) + V(i, 1:end-1));
                n_sub(i,:) = 0.5*(n(i, 2:end) + n(i, 1:end-1));
                p_sub(i,:) = 0.5*(p(i, 2:end) + p(i, 1:end-1));
                c_sub(i,:) = 0.5*(c(i, 2:end) + c(i, 1:end-1));
                a_sub(i,:) = 0.5*(a(i, 2:end) + a(i, 1:end-1));
                
                dVdx(i,:) = (V(i, 2:end) - V(i, 1:end-1))./(x(2:end) - x(1:end-1));   
                dndx(i,:) = (n(i, 2:end) - n(i, 1:end-1))./(x(2:end) - x(1:end-1));
                dpdx(i,:) = (p(i, 2:end) - p(i, 1:end-1))./(x(2:end) - x(1:end-1));
                dcdx(i,:) = (c(i, 2:end) - c(i, 1:end-1))./(x(2:end) - x(1:end-1));
                dadx(i,:) = (a(i, 2:end) - a(i, 1:end-1))./(x(2:end) - x(1:end-1));              
            end
                
            % Diffusion coefficients
            if par.prob_distro_function == 'Fermi'
                for jj = 1:length(x)
                    Dn(i,jj) = distro_fun.D(n(i,jj), dev.Dnfun(jj,:), dev.n_fd(jj,:));
                    Dp(i,jj) = distro_fun.D(p(i,jj), dev.Dpfun(jj,:), dev.p_fd(jj,:));
                end
            end
        
            if par.prob_distro_function == 'Boltz'
                Dn_mat = mue_mat*par.kB*par.T;
                Dp_mat = muh_mat*par.kB*par.T;
            end

            % Particle fluxes (remember F = -dVdx)
            jdd.ndrift = mue_mat.*n_sub.*(dVdx - gradEA_mat);
            jdd.ndiff = -Dn_mat.*(dndx - ((n_sub./Nc_mat).*gradNc_mat));
            jdd.pdrift = muh_mat.*p_sub.*(-dVdx + gradIP_mat);
            jdd.pdiff = -Dp_mat.*(dpdx - ((p_sub./Nv_mat).*gradNv_mat));
                     
            switch par.N_ionic_species
                case 0
                    jdd.cdrift = zeros(length(t), length(xout));
                    jdd.cdiff = zeros(length(t), length(xout));
                    jdd.adrift = zeros(length(t), length(xout));
                    jdd.adiff = zeros(length(t), length(xout));
                case 1
                    jdd.cdrift = mu_cat.*c_sub.*-dVdx;
                    jdd.cdiff = -mu_cat.*par.kB*par.T.*dcdx;
                    jdd.adrift = zeros(length(t), length(xout));
                    jdd.adiff = zeros(length(t), length(xout));               
                case 2
                    jdd.cdrift = mu_cat.*c_sub.*-dVdx;
                    jdd.cdiff = -mu_cat.*par.kB*par.T.*dcdx;
                    jdd.adrift = -mu_ani.*a_sub.*-dVdx;
                    jdd.adiff = -mu_ani.*par.kB*par.T.*dadx;
            end
            
            % Note these have a negative sign compared with the expressions in DFPDE
            jdd.n = par.mobset*(jdd.ndrift + jdd.ndiff );
            jdd.p = par.mobset*(jdd.pdrift + jdd.pdiff);
            jdd.a = par.mobseti*(jdd.adrift + jdd.adiff);
            jdd.c = par.mobseti*(jdd.cdrift + jdd.cdiff);

            % Displacement current
            [~, dFdt] = gradient(-dVdx, xout, t);
            j.disp = par.epp0.*eppmat.*dFdt;
            
            jdd.disp = j.disp;
            % The total flux here includes the sign of the carrier
            jdd.tot = -jdd.n + jdd.p - jdd.a + jdd.c + jdd.disp;

            Jdd.ndrift = jdd.ndrift*par.e;
            Jdd.ndiff = jdd.ndiff*par.e;
            Jdd.pdrift = jdd.pdrift*par.e;
            Jdd.pdiff = jdd.pdiff*par.e;
            Jdd.adrift = jdd.adrift*par.e;
            Jdd.adiff = jdd.adiff*par.e;
            Jdd.cdrift = jdd.cdrift*par.e;
            Jdd.cdiff = jdd.cdiff*par.e;

            Jdd.n = par.e*jdd.n;
            Jdd.p = par.e*jdd.p;
            Jdd.a = par.e*jdd.a;
            Jdd.c = par.e*jdd.c;
            Jdd.disp = par.e*jdd.disp;
            Jdd.tot = par.e*jdd.tot;
        end

        function [FV, Frho] = calcF(sol)
            % Electric field caculation
            % FV = Field calculated from the gradient of the potential
            % Frho = Field calculated from integrated space charge density
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            FV = -gradient(V, x, t);                      % Electric field calculated from V

            if nargout > 1
                rho = dfana.calcrho(sol);
                Frho = cumtrapz(x, rho, 2)./(dev.epp.*par.epp0) + FV(:,1);
            end
        end

        function [FV, Frho] = calcF_ihalf(sol)
            % Electric field caculation
            % FV = Field calculated from the gradient of the potential
            % Frho = Field calculated from integrated space charge density
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            for i=1:length(t)
                [~, dVdx_ihalf(i,:)] = pdeval(0,x,V(i,:),par.x_ihalf);
            end
            FV = -dVdx_ihalf;

            if nargout > 1
                rho = dfana.calcrho_ihalf(sol);
                Frho = cumtrapz(par.x_ihalf, rho, 2)./(par.dev_ihalf.epp.*par.epp0) + FV(:,1);
            end
        end

        function rho = calcrho(sol)
            % Calculates the space charge density
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            % charge density
            rho = -n + p - a + c - dev.NA + dev.ND + dev.Nani - dev.Ncat;
        end

        function rho = calcrho_ihalf(sol)
            % Calculates the space charge density
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            n_ihalf = getvarihalf(n);
            p_ihalf = getvarihalf(p);
            a_ihalf = getvarihalf(a);
            c_ihalf = getvarihalf(c);

            % charge density
            rho = -n_ihalf + p_ihalf - a_ihalf + c_ihalf - par.dev_ihalf.NA + par.dev_ihalf.ND + par.dev_ihalf.Nani - par.dev_ihalf.Ncat;
        end

        function Vapp = calcVapp(sol)
            [~,t,~,par,~,~,~,~,~,~] = dfana.splitsol(sol);
            switch par.V_fun_type
                case 'constant'
                    Vapp = ones(1,length(t))*par.Vapp;
                otherwise
                    Vapp_fun = fun_gen(par.V_fun_type);
                    Vapp = Vapp_fun(par.V_fun_arg, t);
            end
        end

        function stats = JVstats(JVsol)
            % A function to pull statistics from a JV sweep using DOJV
            % JVsol - a solution from DOJV
            if isfield(JVsol, 'ill')
                if isfield(JVsol.ill, 'f')
                    Vapp = dfana.calcVapp(JVsol.ill.f);
                    Vapp = Vapp';
                    J = dfana.calcJ(JVsol.ill.f);
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
                    Vapp = Vapp';
                    J = dfana.calcJ(JVsol.ill.r);
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

        function value = PLt(sol)
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            Bmat = dev.B;
            value = trapz(x,(dev.B.*(n.*p-dev.ni.^2)),2);
        end

        function VQFL = calcVQFL(sol)
            % Get QFLs
            [~, ~, Efn, Efp] = dfana.QFLs(sol);
            par = sol.par;
            if par.pleft >= par.nleft && par.nright >= par.pright
                % p-type left boundary, n-type right boundary
                VQFL = Efn(:, end) - Efp(:, 1);
            elseif par.nleft >= par.nright && par.pright >= par.nright
                % n-type left boundary, p-type right boundary
                VQFL = Efn(:, 1) - Efp(:, end);
            else
                % If all equal the choose arbitrary boundaries
                VQFL = Efn(:, end) - Efp(:, 1);
            end
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

        function sigma_ion = calcsigma_ion(sol)
            % calculates the integrated space charge density
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            rho_ion = c-a;
            sigma_ion = trapz(x, rho_ion, 2);
        end

        function Fion = calcFion(sol)
           [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

           rhoion = c - dev.Ncat - a + dev.Nani;
           Fion = cumtrapz(x, rhoion, 2)./(dev.epp*par.epp0);

        end

        function Vion = calcVion(sol)
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            Fion = dfana.calcFion(sol);
            Vion = -cumtrapz(x, Fion,2);
        end

        function [U,Ux] = pdentrp(singular,m,xL,uL,xR,uR,xout)
            % PDENTRP  Interpolation helper function for PDEPE.
            %   [U,UX] = PDENTRP(M,XL,UL,XR,UR,XOUT) uses solution values UL at XL and UR at XR
            %   for successive mesh points XL < XR to interpolate the solution values U and
            %   the partial derivative with respect to x, UX, at arguments XOUT(i) with
            %   XL <= XOUT(i) <= XR.  UL and UR are column vectors. Column i of the output
            %   arrays U, UX correspond to XOUT(i).
            %
            %   See also PDEPE, PDEVAL, PDEODES.

            %   Lawrence F. Shampine and Jacek Kierzenka
            %   Copyright 1984-2013 The MathWorks, Inc.
            %   $Revision: 1.5.4.4.54.1 $  $Date: 2013/09/27 03:10:22 $

            xout = xout(:)';
            nout = length(xout);

            U = uL(:,ones(1,nout));
            Ux = zeros(size(U));

            uRL = uR - uL;
            % Use singular interpolant on all subintervals.
            if singular
                U  = U + uRL*((xout .^ 2 - xL^2) / (xR^2 - xL^2));
                Ux =     uRL*(2*xout / (xR^2 - xL^2));
            else
                switch m
                    case 0
                        U  = U + uRL*( (xout - xL) / (xR - xL));
                        Ux =     uRL*(ones(1,nout) / (xR - xL));
                    case 1
                        U  = U + uRL*(log(xout/xL) / log(xR/xL));
                        Ux =     uRL*( (1 ./ xout) / log(xR/xL));
                    case 2
                        U  = U + uRL*((xR ./ xout) .* ((xout - xL)/(xR - xL)));
                        Ux =     uRL*((xR ./ xout) .* (xL ./ xout)/(xR - xL));
                end
            end
        end
    end
end
