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

            switch par.N_ionic_species
                case 0
                    c = zeros(length(t), length(x));
                    a = zeros(length(t), length(x));
                    dev.Ncat = zeros(1, length(x));
                    dev.Nani = zeros(1, length(x));
                    par.dev.Ncat = zeros(1, length(x));
                    par.dev.Nani = zeros(1, length(x));
                    par.dev_sub.Ncat = zeros(1, length(x) - 1);
                    par.dev_sub.Nani = zeros(1, length(x) - 1);
                case 1
                    c = u(:,:,4);
                    a = zeros(length(t), length(x));
                    dev.Nani = zeros(1, length(x));
                    par.dev.Nani = zeros(1, length(x));
                    par.dev_sub.Nani = zeros(1, length(x) - 1);
                case 2
                    c = u(:,:,4);
                    a = u(:,:,5);
            end
        end

        function [Ecb, Evb, Efn, Efp] = calcEnergies(sol)
            % u is the solution structure
            % Simple structure names
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            Ecb = dev.Phi_EA-V;                                 % Conduction band potential
            Evb = dev.Phi_IP-V;                                 % Valence band potential
            Efn = zeros(size(n,1), size(n,2));
            Efp = zeros(size(n,1), size(n,2));

            switch par.prob_distro_function
                case 'Fermi'

                    for i = 1:size(n,1)           % time
                        for j = 1:size(n,2)       % position
                            Efn(i,j) = distro_fun.Efn_fd_fun(n(i,j), dev.Efn(j,:),  dev.n_fd(j,:));
                            Efp(i,j) = distro_fun.Efp_fd_fun(p(i,j), dev.Efp(j,:),  dev.p_fd(j,:));
                        end
                    end
                    Efn = Efn-V;
                    Efp = Efp-V;

                case 'Blakemore'
                    Efn = real(Ecb + (par.kB*par.T/par.q)*log(n./(dev.Nc - par.gamma*n)));
                    Efp = real(Evb - (par.kB*par.T/par.q)*log(p./(dev.Nv - par.gamma*p)));

                case 'Boltz'
                    Efn = real(Ecb + (par.kB*par.T/par.q)*log(n./dev.Nc));        % Electron quasi-Fermi level
                    Efp = real(Evb - (par.kB*par.T/par.q)*log(p./dev.Nv));        % Hole quasi-Fermi level
            end

        end

        function [J, j, x] = calcJ(sol)
            % Current, J and flux, j calculation from continuity equations
            % obtain SOL components for easy referencing
            [u,t,xmesh,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            n_sub = getvar_sub(n);
            p_sub = getvar_sub(p);
            a_sub = getvar_sub(a);
            c_sub = getvar_sub(c);

            x = par.x_sub;
            [~,~,g] = dfana.calcg(sol);

            [~, dndt] = gradient(n_sub, x, t);
            [~, dpdt] = gradient(p_sub, x, t);
            [~, dadt] = gradient(a_sub, x, t);
            [~, dcdt] = gradient(c_sub, x, t);

            % Recombination
            r = dfana.calcr(sol, "sub");

            djndx = -dndt + g - r.tot;
            djpdx = -dpdt + g - r.tot;
            djadx = -dadt;                  % Add source terms as necessary
            djcdx = -dcdt;                  % Add source terms as necessary

            deltajn = cumtrapz(x, djndx, 2);
            deltajp = cumtrapz(x, djpdx, 2);
            deltaja = cumtrapz(x, djadx, 2);
            deltajc = cumtrapz(x, djcdx, 2);

            %% Currents from the boundaries
            jn_l = -par.sn_l*(n(:, 1) - par.n0_l);
            jn_r = par.sn_r*(n(:, end) - par.n0_r);

            jp_l = -par.sp_l*(p(:, 1) - par.p0_l);
            jp_r = par.sp_r*(p(:, end) - par.p0_r);

            jc_l = 0;
            jc_r = 0;

            ja_l = 0;
            ja_r = 0;

            % Calculate total electron and hole currents from fluxes
            % Use the minority carrier flux as the boundary condition
            if par.p0_l == par.n0_l && par.n0_r == par.p0_r
                % Intrinsic both sides then integrate from side with lowest
                % density
                if par.n0_r > par.n0_l
                    j.n = jn_l + deltajn;
                else
                    j.n = jn_r + (deltajn-deltajn(:,end));
                end

                if par.p0_r > par.p0_l
                    j.p = jp_l + deltajp;
                else
                    j.p = jp_r + (deltajp - deltajp(:,end));
                end
            elseif par.p0_l >= par.n0_l && par.n0_r >= par.p0_r
                % p-type left boundary, n-type right boundary
                j.n = jn_l + deltajn;
                j.p = jp_r + (deltajp - deltajp(:,end));
            elseif par.n0_l >= par.n0_r && par.p0_r >= par.n0_r
                % n-type left boundary, p-type right boundary
                j.n = jn_r + (deltajn-deltajn(:,end));
                j.p = jp_l + deltajp;
            elseif par.p0_l >= par.n0_l && par.p0_r >= par.n0_r...
                    || par.n0_l >= par.p0_l && par.n0_r >= par.p0_r
                % p-type both boundaries or n-type both boundaries
                j.n = jn_l + deltajn;
                j.p = jp_l + deltajp;
            end

            j.c = jc_l + deltajc;
            j.a = ja_l + deltaja;

            % Apply switches and accelerators
            j.n = par.mobset*j.n;
            j.p = par.mobset*j.p;
            j.c = par.mobseti*par.K_c*j.c;
            j.a = par.mobseti*par.K_a*j.a;

            % displacement flux
            FV_sub = dfana.calcF(sol, "sub");

            [~, FV_sub_dt] = gradient(FV_sub, x, t);
            j.disp = par.epp0.*par.dev_sub.epp.*FV_sub_dt;

            J.n = j.n*-par.e;
            J.p = j.p*par.e;
            J.c = j.c*par.z_c*par.e;
            J.a = j.a*par.z_a*par.e;
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
        
        function Pin = calcPin(sol)
            % Incident optical power density
            % Note this integrates across the available spectrum 
            % within AM15.xls 
            if strcmp(sol.par.optical_model, 'Beer-Lambert')
                AM15_data = xlsread('AM15.xls');
                Pin = 1e-3*trapz(AM15_data(:,1), AM15_data(:,2));
            else
                warning('No incident photon spectrum available, assuming Pin = 0.1 W cm-2')
                Pin = 0.1;
            end
        end

        
        function [r, ns, ps, alpha_xn, beta_xp] = calcr(sol, mesh_option)
            % Calculate the recombination rate on i-half mesh
            % obtain SOL components for easy referencing
            % MESH_OPTION = "whole" for input mesh or "sub" for subinetrval
            % mesh
            [u,t,x_input,par,~,n,p,a,c,V] = dfana.splitsol(sol);

            switch mesh_option
                case "whole"
                    dev = par.dev;
                    x = x_input;
                    n = u(:,:,2);
                    p = u(:,:,3);
                case "sub"
                    dev = par.dev_sub;
                    x = par.x_sub;
                    n = getvar_sub(u(:,:,2));
                    p = getvar_sub(u(:,:,3));
            end

            dVdx = zeros(length(t), length(x));
            for i=1:length(t)
                [~, dVdx(i,:)] = pdeval(0, x_input, V(i,:), x);
            end
            vsr_zone = repmat(dev.vsr_zone, length(t), 1);
            srh_zone = repmat(dev.srh_zone, length(t), 1);

            xprime_n = dev.xprime_n;
            xprime_p = dev.xprime_p;
            sign_xn = repmat(dev.sign_xn, length(t), 1);    % 1 if xn increasing, -1 if decreasing wrt x
            sign_xp = repmat(dev.sign_xp, length(t), 1);    % 1 if xp increasing, -1 if decreasing wrt x
            alpha0_xn = repmat(dev.alpha0_xn, length(t), 1);
            beta0_xp = repmat(dev.beta0_xp, length(t), 1);

            alpha_xn = (sign_xn.*par.q.*dVdx./(par.kB*par.T)) + alpha0_xn;
            beta_xp = (sign_xp.*par.q.*-dVdx./(par.kB*par.T)) + beta0_xp;

            % Band-to-band
            r.btb = dev.B.*(n.*p - dev.ni.^2);
            % Bulk SRH
            r.srh = srh_zone.*(n.*p - dev.ni.^2)...
                ./(dev.taun.*(p+dev.pt) + dev.taup.*(n+dev.nt));
            % Volumetric surface SRH
            ns = n.*exp(-alpha_xn.*xprime_n); % Projected electron surface density
            ps = p.*exp(-beta_xp.*xprime_p);  % Projected hole surface density
            r.vsr = vsr_zone.*(ns.*ps - dev.ni.^2)...
                ./(dev.taun_vsr.*(ps + dev.pt) + dev.taup_vsr.*(ns + dev.nt));
            % Total
            r.tot = r.btb + r.srh + r.vsr;
        end

        function [Jdd, jdd, xout] = calcJdd(sol)
            % Calculates drift and diffusion currents at every point and all times -
            % NOTE: UNRELIABLE FOR TOTAL CURRENT as errors in the calculation of the
            % spatial gradients mean that the currents do not cancel properly
            % obtain SOL components for easy referencing
            [~,t,x,par,~,n,p,a,c,V] = dfana.splitsol(sol);
            xout = par.x_sub;
            dev = par.dev_sub;

            % Property matrices
            eppmat = dev.epp;
            mu_n_mat = dev.mu_n;
            mu_p_mat = dev.mu_p;
            mu_cat = dev.mu_c;
            mu_ani = dev.mu_a;
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
            switch par.prob_distro_function
                case 'Fermi'
                    for jj = 1:length(x)
                        Dn_mat(i,jj) = distro_fun.D(n(i,jj), dev.Dnfun(jj,:), dev.n_fd(jj,:));
                        Dp_mat(i,jj) = distro_fun.D(p(i,jj), dev.Dpfun(jj,:), dev.p_fd(jj,:));
                    end

                case 'Blakemore'
                    Dn_mat = mu_n_mat.*par.kB.*par.T.*(Nc_mat./(Nc_mat - par.gamma.*n_sub));
                    Dp_mat = mu_p_mat.*par.kB.*par.T.*(Nv_mat./(Nv_mat - par.gamma.*p_sub));

                case 'Boltz'
                    Dn_mat = mu_n_mat*par.kB*par.T;
                    Dp_mat = mu_p_mat*par.kB*par.T;
            end

            % Particle fluxes (remember F = -dVdx)
            jdd.ndrift = mu_n_mat.*n_sub.*(dVdx - gradEA_mat);
            jdd.ndiff = -Dn_mat.*(dndx - ((n_sub./Nc_mat).*gradNc_mat));
            jdd.pdrift = mu_p_mat.*p_sub.*(-dVdx + gradIP_mat);
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

        function [FV, Frho] = calcF(sol, mesh_option)
            % Electric field caculation
            % FV = Field calculated from the gradient of the potential
            % Frho = Field calculated from integrated space charge density
            [u,t,x_whole,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            switch mesh_option
                case "whole"
                    x = x_whole;
                case "sub"
                    x = par.x_sub;
            end

            for i=1:length(t)
                [~, dVdx(i,:)] = pdeval(0,x_whole,V(i,:),x);
            end
            FV = -dVdx;

            if nargout > 1
                rho = dfana.calcrho(sol, "sub");
                Frho = cumtrapz(x, rho, 2)./(par.dev_sub.epp.*par.epp0) + FV(:,1);
            end
        end

        function rho = calcrho(sol, mesh_option)
            % Calculates the space charge density
            [u,t,x,par,dev_in,n_whole,p_whole,a_whole,c_whole,V_whole] = dfana.splitsol(sol);

            switch mesh_option
                case "whole"
                    dev = par.dev;
                    n = n_whole;
                    p = p_whole;
                    a = a_whole;
                    c = c_whole;
                case "sub"
                    dev = par.dev_sub;
                    n = getvar_sub(n_whole);
                    p = getvar_sub(p_whole);
                    a = getvar_sub(a_whole);
                    c = getvar_sub(c_whole);
            end

            NA = repmat(dev.NA, length(t), 1);
            ND = repmat(dev.ND, length(t), 1);
            Nani = repmat(dev.Nani, length(t), 1);
            Ncat = repmat(dev.Ncat, length(t), 1);
            % charge density
            rho = -n + p - NA + ND + par.z_a*a + par.z_c*c - par.z_c*Nani - par.z_c*Ncat;
        end

        function Vapp = calcVapp(sol)
            [~,t,~,par,~,~,~,~,~,~] = dfana.splitsol(sol);
            switch par.V_fun_type
                case 'constant'
                    Vapp = ones(1,length(t))*par.V_fun_arg(1);
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

        function value = calcPLt(sol)
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

            Bmat = dev.B;
            value = trapz(x,(dev.B.*(n.*p-dev.ni.^2)),2);
        end

        function DeltaQFL = calcDeltaQFL(sol)
            % Get QFLs
            [~, ~, Efn, Efp] = dfana.calcEnergies(sol);
            par = sol.par;
            if par.p0_l >= par.n0_l && par.n0_r >= par.p0_r
                % p-type left boundary, n-type right boundary
                DeltaQFL = Efn(:, end) - Efp(:, 1);
            elseif par.n0_l >= par.n0_r && par.p0_r >= par.n0_r
                % n-type left boundary, p-type right boundary
                DeltaQFL = Efn(:, 1) - Efp(:, end);
            else
                % If all equal the choose arbitrary boundaries
                DeltaQFL = Efn(:, end) - Efp(:, 1);
            end
        end

        function deltaV = deltaVt(sol, p1, p2)
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            % Calculates the electrostatic potential difference as a function of time
            % between two points P1 and P2
            deltaV = V(:,p1) - V(:,p2);
        end

        function sigma = calcsigma(sol)
            % calculates the integrated space charge density
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            rho = dfana.calcrho(sol, "whole");
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
