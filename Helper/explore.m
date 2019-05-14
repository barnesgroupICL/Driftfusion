classdef explore

    methods (Static)

        function parexsol = explore2par(par_base, parnames, parvalues, JVpnts)
            % EXPLORE2PAR is used to explore 2 different parameters using the parallel computing toolbox.
            % The code is likely to require modification for individual parameters
            % owing to possible dependencies.
            % PAR_BASE is the base parameter set
            % PARNAMES is a cell array with the parameter names in - check these
            % carefully to avoid heartache later
            % PARVALUES is matrix with the parameter value ranges e.g.
            tic
            disp('Starting parameter exploration');
            disp(['Parameter 1: ', parnames(1)]);
            disp(['Parameter 2: ', parnames(2)]);
            parval1 = cell2mat(parvalues(1));
            parval2 = cell2mat(parvalues(2));
            str1 = char(parnames(1));
            str2 = char(parnames(2));

            j = 1;
            parfor i = 1:length(parval1)
                
                try
                par = par_base;
                par.Ana = 0;
                par = explore.helper(par, str1, parval1(i));
                
                % If the parameter is 
                if strmatch('dcell', parnames(1)) ~= 0
                    % sets PCELL at the required position to provide a point
                    % density of one point per nanometer
                    layerpoints = round(parval1(i)*1e7);
                    par = explore.helper(par, ['p', str1(2:end)], layerpoints);
                    par.pcell(1,4) = layerpoints*1;
                end

                % Rebuild device
                par.xx = pc.xmeshini(par);
                par.dev = pc.builddev(par);

                Voc_f = zeros(1, length(parval2));
                Voc_r = zeros(1, length(parval2));
                Jsc_f = zeros(1, length(parval2));
                Jsc_r = zeros(1, length(parval2));
                mpp_f = zeros(1, length(parval2));
                mpp_r = zeros(1, length(parval2));
                FF_f = zeros(1, length(parval2));
                FF_r = zeros(1, length(parval2));
                n_av = zeros(1, length(parval2));
                p_av = zeros(1, length(parval2));
                n_f = zeros(length(parval2), 2000);   % final electron profile of stabilised sol
                p_f = zeros(length(parval2), 2000);
                a_f = zeros(length(parval2), 2000);
                V_f = zeros(length(parval2), 2000);
                x = zeros(length(parval2), 2000);
                Voc_stable = zeros(length(parval2), JVpnts);
                PLint = zeros(length(parval2), JVpnts);
                Vapp_f = zeros(1, JVpnts);
                J_f = zeros(length(parval2), JVpnts);
                Vapp_r = zeros(1, JVpnts);
                J_r = zeros(length(parval2), JVpnts);

                soleq = equilibrate(par);

                for j = 1:length(parval2)

                    runN = (i-1)*length(parval2) + j;
                    disp(['Run no. ', num2str(runN), ', ', str1 ,' = ', num2str(parval1(i)), ' , ',str2,' = ', num2str(parval2(j))]);
                    
                    % If the second parameter name is intensity then set
                    % INT using PARVAL2 else set the appropriate
                    % parameter to the required value in par and Int is
                    % assumed to be 1
                    if strmatch('Int', parnames(2)) == 1
                        Int = parval2(j);
                    else
                        par = explore.helper(par, str2, parval2(j));
                        Int = 1;
                        % Rebuild device
                        par.xx = pc.xmeshini(par);
                        par.dev = pc.builddev(par);
                    end

                    JV = doJV(soleq.ion, 50e-3, JVpnts, Int, 1, 0, 1.3, 2);
                    stats = dfana.JVstats(JV);

                    Vapp_f(j,:) = JV.ill.f.Vapp;
                    J_f(j,:) =  JV.ill.f.J.tot(:,end);
                    Vapp_r(j,:) = JV.ill.r.Vapp;
                    J_r(j,:) = JV.ill.r.J.tot(:,end);
                    Voc_f(j) = stats.Voc_f;
                    Voc_r(j) = stats.Voc_r;
                    Jsc_f(j) = stats.Jsc_f;
                    Jsc_r(j) = stats.Jsc_r;
                    mpp_f(j) = stats.mpp_f;
                    mpp_r(j) = stats.mpp_r;
                    FF_f(j) = stats.FF_f;
                    FF_r(j) = stats.FF_r;
                    
                    % For steady-state solution
                    sol_ill = lighton_Rs(soleq.ion, Int, 1, 1e-6, 0, JVpnts);

                    Voc_stable(j,:) = dfana.Voct(sol_ill);
                    PLint(j,:) = dfana.PLt(sol_ill);
                    n_av(j) = mean(sol_ill.u(end, par.pcum(2):par.pcum(5),1));
                    p_av(j) = mean(sol_ill.u(end, par.pcum(2):par.pcum(5),2));
                    n_f = explore.writevar(n_f, j, par.xx, sol_ill.u(end,:,1));
                    p_f = explore.writevar(p_f, j, par.xx, sol_ill.u(end,:,2));
                    a_f = explore.writevar(a_f, j, par.xx, sol_ill.u(end,:,3));
                    V_f = explore.writevar(V_f, j, par.xx, sol_ill.u(end,:,4));
                    x = explore.writevar(x, j, par.xx, sol_ill.x);
                end

                A(i,:) = Voc_f;
                B(i,:) = Voc_r;
                C(i,:) = Jsc_f;
                D(i,:) = Jsc_r;
                E(i,:) = mpp_f;
                F(i,:) = mpp_r;
                G(i,:) = FF_f;
                H(i,:) = FF_r;
                J(i,:,:) = Voc_stable;
                K(i,:,:) = PLint;
                AA(i,:,:) = Vapp_f;
                BB(i,:,:) = J_f;
                CC(i,:,:) = Vapp_r;
                DD(i,:,:) = J_r;
                EE(i,:) = n_av;
                FF(i,:) = p_av;
                GG(i,:,:) = n_f;
                HH(i,:,:) = p_f;
                II(i,:,:) = a_f;
                JJ(i,:,:) = V_f;
                KK(i,:,:) = x;
                
                catch
                   warning(['Run no. ', num2str(runN), ', ', str1 ,' = ', num2str(parval1(i)), ' , ',str2,' = ', num2str(parval2(j)), ' failed']);
                    
                end
            end

            % Store solutions in output struct
            parexsol.stats.Voc_f = A;
            parexsol.stats.Voc_r = B;
            parexsol.stats.Jsc_f = C;
            parexsol.stats.Jsc_r = D;
            parexsol.stats.mpp_f = E;
            parexsol.stats.mpp_r = F;
            parexsol.stats.FF_f = G;
            parexsol.stats.FF_r = H;
            parexsol.Vapp_f = AA;
            parexsol.J_f = BB;
            parexsol.Vapp_r = CC;
            parexsol.J_r = DD;
            parexsol.stats.Voc_stable = J;
            parexsol.stats.PLint = K;
            parexsol.parnames = parnames;
            parexsol.parvalues = parvalues;
            parexsol.parval1 = parval1;
            parexsol.parval2 = parval2;
            parexsol.par_base = par_base;
            parexsol.stats.n_av = EE;
            parexsol.stats.p_av = FF;
            parexsol.n_f = GG;
            parexsol.p_f = HH;
            parexsol.a_f = II;
            parexsol.V_f = JJ;
            parexsol.x = KK;

            toc

        end

        function var = writevar(var, j, xx, arr)
            % EXPLORE.WRITEVAR writes array ARR to variable VAR using variable length
            % assignment, which is not usually allowed in PARFOR loops. This
            % allows the solution with different length vectors to be stored.
            var(j,1:length(xx)) = arr;
        end

        function par = helper(par, parname, parvalue)
            % takes parameter set and sets parname to parvalue- workaround for parallel
            % computing loops
            eval(['par.',parname,'=parvalue']);
        end

        function plotPL(parexsol)

            figure(3000)
            s1 = surf(parexsol.parval1, parexsol.parval2, parexsol.stats.PLint);
            ylabel(parexsol.parnames(1))
            xlabel(parexsol.parnames(2))
            set(s1,'YScale','log');
            zlabel('PL intensity [cm-2s-1]')
            shading interp
            colorbar

        end

        function plotVoc(parexsol, xlogon, ylogon, zlogon)

            figure(3001)
            surf(parexsol.parval2, parexsol.parval1+60e-7, parexsol.stats.Voc_f)
            s1 = gca;
            %ylabel('Ion density [cm-3]')
            ylabel(parexsol.parnames{1,1})
            xlabel(parexsol.parnames{1,2})
            zlabel('Voc F scan [V]')
            xlim([parexsol.parval2(1), parexsol.parval2(end)]);
            ylim([parexsol.parval1(1)+60e-7, parexsol.parval1(end)+60e-7])

            %caxis([0.75, 0.95])
            if xlogon
                set(s1,'XScale','log');
            else
                set(s1,'XScale','linear');
            end
            if ylogon
                set(s1,'YScale','log');
            else
                set(s1,'YScale','linear');
            end
            zlabel('taurise [s]')
            shading interp
            colorbar
            cb = colorbar();
            if zlogon
                cb.Ruler.Scale = 'log';
                cb.Ruler.MinorTick = 'on';
            end

            figure(3002)
            surf(parexsol.parval2, parexsol.parval1+60e-7, parexsol.stats.Voc_r)
            s1 = gca;
            %ylabel('Ion density [cm-3]')
            ylabel(parexsol.parnames{1,1})
            xlabel(parexsol.parnames{1,2})
            zlabel('Voc R scan [V]')
            xlim([parexsol.parval2(1), parexsol.parval2(end)]);
            ylim([parexsol.parval1(1)+60e-7, parexsol.parval1(end)+60e-7])
            set(s1,'YScale','log');
            shading interp
            colorbar
            %caxis([0.75, 0.95])
            if xlogon
                set(s1,'XScale','log');
            else
                set(s1,'XScale','linear');
            end
            if ylogon
                set(s1,'YScale','log');
            else
                set(s1,'YScale','linear');
            end
            zlabel('taurise [s]')
            shading interp
            colorbar
            cb = colorbar();
            if zlogon
                cb.Ruler.Scale = 'log';
                cb.Ruler.MinorTick = 'on';
            end

        end

        function plotstat_2D_parval1(parexsol, yproperty, logx, logy)
            % YPROPERTY is string with the property name. Properties must
            % be one of those contained in the PAREXSOL.STATS structure
            eval(['y = parexsol.stats.', yproperty])
            for i=1:length(parexsol.parval1)
                figure(3010)
                if logx == 0 && logy == 0
                    plot(parexsol.parval2, y(i,:))
                elseif logx == 1 && logy == 0
                    semilogx(parexsol.parval2, y(i,:))
                elseif logx == 0 && logy == 1
                    semilogy(parexsol.parval2, y(i,:))
                elseif logx == 1 && logy == 1
                    loglog(parexsol.parval2, y(i,:))
                end
                hold on
            end
            xlabel(parexsol.parnames{1,2})
            ylabel(yproperty)
            legstr = (num2str((parexsol.parval1'+60e-7)*1e7));
            legend(legstr)
            %xlim([min(parexsol.parval2), max(parexsol.parval2)])
            xlim([1e-3, 10])
            hold off
        end

        function plotstat_2D_parval2(parexsol, yproperty, logx, logy)
            % YPROPERTY is string with the property name. Properties must
            % be one of those contained in the PAREXSOL.STATS structure

            figure(3010)

            eval(['y = parexsol.stats.', yproperty,';'])

            if strmatch('Voc_stable', yproperty) ==1
               y = y(:,:,end);
            end

            if strmatch('Jsc_r', yproperty) ==1
               y = -y;
            end

            x = (parexsol.parval1'+ 60e-7)*1e7;

            for i=1:length(parexsol.parval2)
                            % Rebuild solutions
                            if logx == 0 && logy == 0
                                plot(x, y(:,i))
                            elseif logx == 1 && logy == 0
                                semilogx(x, y(:,i))
                            elseif logx == 0 && logy == 1
                                semilogy(x, y(:,i))
                            elseif logx == 1 && logy == 1
                                loglog(x, y(:,i))
                            end
                            hold on

            end
            %xlabel(parexsol.parnames{1,1})
            xlabel('Active layer thickness [nm]')
            ylabel(yproperty)
            legstr = (num2str(parexsol.parval2'));
            legend(legstr)
            xlim([min(x), max(x)])

            hold off
        end

        function plotEL(parexsol)
            % plots final energy level diagrams
            figure(3011)
            for i=1:length(parval1)
                for j = 1:length(parval2)
                    % Rebuild solutions
                    sol.u(1,:,1) = parexsol.nf(i, j, :);
                    sol.u(1,:,2) = parexsol.pf(i, j, :);
                    sol.u(1,:,3) = parexsol.af(i, j, :);
                    sol.u(1,:,4) = parexsol.Vf(i, j, :);
                    sol.t = 0;
                    sol.x = parexsol.x(i,j,:);
                    sol.par = parexsol.par_base;

                end
            end
        end

        function plotprof_2D(parexsol, yproperty, par1logical, par2logical,logx,logy)
            eval(['y = parexsol.', yproperty,';']);
            if length(par1logical) > length(parexsol.parval1)
                par1logical = par1logical(1:length(parexsol.parval1));
            end

            if length(par2logical) > length(parexsol.parval2)
                par2logical = par2logical(1:length(parexsol.parval2));
            end

            % Replace zeros with NaNs
            y(y==0) = NaN;
            parexsol.x(parexsol.x==0) = NaN;

            if strmatch(yproperty,'a_f')
                y = y-parexsol.par_base.Nion(1);
            end

            figure(3012)
            for i=1:length(parexsol.parval1)
                if par1logical(i) == 1
                    for j = 1:length(parexsol.parval2)
                        if par2logical(j) == 1
                            % Rebuild solutions
                            if logx
                                semilogx(squeeze(parexsol.x(i, j, :).*1e7), squeeze(y(i, j, :)));
                            elseif logy
                                semilogy(squeeze(parexsol.x(i, j, :).*1e7), squeeze(y(i, j, :)));
                            elseif logx ==1 && logy ==1
                                loglog(squeeze(parexsol.x(i, j, :).*1e7), squeeze(y(i, j, :)));
                            else
                                plot(squeeze(parexsol.x(i, j, :).*1e7), squeeze(y(i, j, :)));
                            end
                            hold on
                        end
                    end
                end
            end
            hold off
            xlabel('Position [nm]')
            ylabel(yproperty)

        end

        function plotU(parexsol, par1logical, par2logical,logx,logy)

            if length(par1logical) > length(parexsol.parval1)
                par1logical = par1logical(1:length(parexsol.parval1));
            end

            if length(par2logical) > length(parexsol.parval2)
                par2logical = par2logical(1:length(parexsol.parval2));
            end

            UHTL = zeros(length(par1logical),length(par2logical));
            UHTLint = zeros(length(par1logical),length(par2logical));
            Ubulk = zeros(length(par1logical),length(par2logical));
            UETLint = zeros(length(par1logical),length(par2logical));
            UETL = zeros(length(par1logical),length(par2logical));
            Utotal = zeros(length(par1logical),length(par2logical));

            figure(3013)
            par = parexsol.par_base;

            for i=1:length(parexsol.parval1)
                if par1logical(i) == 1
                    for j = 1:length(parexsol.parval2)
                        if par2logical(j) == 1
                            % Rebuild solutions
                            par = explore.helper(par, parexsol.parnames{1,1}, parexsol.parval1(i));

                            if strmatch('dcell(1,4)', parexsol.parnames(1)) ~= 0
                                pcontact = round(parexsol.parval1(i)*1e7);
                                par.pcell(1,4) = pcontact*1;
                            end

                            % Rebuild device
                            par.xx = pc.xmeshini(par);
                            par.dev = pc.builddev(par);

                            dev = par.dev;

                            n = squeeze(parexsol.n_f(i,j,:));
                            p = squeeze(parexsol.p_f(i,j,:));
                            n=n';
                            p=p';
                            n = n(1:length(dev.krad));
                            p = p(1:length(dev.krad));

                            % Recombination
                            Ubtb = dev.krad.*(n.*p - dev.ni.^2);

                            Usrh = ((n.*p - dev.ni.^2)./((dev.taun.*(p+dev.pt)) + (dev.taup.*(n+dev.nt))));

                            U = Ubtb + Usrh;

                            xplot = squeeze(parexsol.x(i, j, :).*1e7);
                            xplot = xplot(1:length(dev.krad));

                            UHTL(i,j) = trapz(xplot(par.pcum0(1):par.pcum0(2)), U(par.pcum0(1):par.pcum0(2)));
                            UHTLint(i,j) = trapz(xplot(par.pcum(1):par.pcum(2)), U(par.pcum(1):par.pcum(2)));                            UHTLint(i,j) = trapz(xplot(par.pcum(1):par.pcum(2)), U(par.pcum(1):par.pcum(2)));
                            Ubulk(i,j) = trapz(xplot(par.pcum(2):par.pcum(5)), U(par.pcum(2):par.pcum(5)));
                            UETLint(i,j) = trapz(xplot(par.pcum(5):par.pcum(6)), U(par.pcum(5):par.pcum(6)));
                            UETL(i,j) = trapz(xplot(par.pcum(6):par.pcum(7)), U(par.pcum(6):par.pcum(7)));
                            Utotal(i,j) = trapz(xplot, U);

                            if logx
                                semilogx(xplot, U);
                            elseif logy
                                semilogy(xplot, U);
                            elseif logx ==1 && logy ==1
                                loglog(xplot, U);
                            else
                                plot(xplot, U);
                            end
                            hold on
                        end
                    end
                end

            end

            figure(3013)
            xlabel('Position [nm]')
            ylabel('U [cm-3s-1]')
            hold off

            for j=1:length(parexsol.parval2)
                if par2logical(j) == 1
                    figure(3014)
                    plot(parexsol.parval1*1e7, UHTL(:,j),...
                        parexsol.parval1*1e7, UHTLint(:,j),...
                        parexsol.parval1*1e7, Ubulk(:,j),...
                        parexsol.parval1*1e7, UETLint(:,j),...
                        parexsol.parval1*1e7, UETL(:,j),...
                        parexsol.parval1*1e7, Utotal(:,j));
                    %hold on
                end
            end
            xlabel('Active layer thickness [nm]')
            ylabel('Recombination rate [cm-3s-1]')
            legend('HTL', 'HTL interface', 'Bulk perovskite', 'ETL interface', 'ETL','Total device')
            hold off

        end


        function plotCE(solex_Voc, solex_eq,logx,logy)
            % plotting function for plotting carrier densities as a
            % function of Voc

            % SOLEX_VOC = EXPLORE open circuit solutions
            % SOLEX_EQ = EXPLORE equilibrium solutions
            % The carrier densities at equilibrium are subtracted from
            % those at Voc to obtain the 'extracted' electronic carrier
            % densities
            % This version requires that the equilibrium solutions are
            % found using a separate run of EXPLORE.EXPLORE2PAR with only a
            % single J element

            par = solex_Voc.par_base;

            for i=1:length(solex_Voc.parval1)   % d_active
                    for j = 1:length(solex_Voc.parval2) % Light intensity
                            % Rebuild solutions
                            par = explore.helper(par, solex_Voc.parnames{1,1}, solex_Voc.parval1(i));

                            if strmatch('dcell(1,4)', solex_Voc.parnames(1)) ~= 0
                                pcontact = round(solex_Voc.parval1(i)*1e7);
                                par.pcell(1,4) = pcontact*1;
                            end

                            % Rebuild device
                            par.xx = pc.xmeshini(par);
                            par.dev = pc.builddev(par);

                            dev = par.dev;

                            n_Voc = squeeze(solex_Voc.n_f(i,j,:))';
                            p_Voc = squeeze(solex_Voc.p_f(i,j,:))';
                            n_eq = squeeze(solex_eq.n_f(i,:));
                            p_eq = squeeze(solex_eq.n_f(i,:));

                            n_CEx(i,:) = n_Voc - n_eq;
                            p_CEx(i,:) = p_Voc - n_eq;

                            n_CE(i,j) = trapz(par.xx(1:length(dev.krad), n_CEx(1:length(dev.krad))));
                            p_CE(i,j) = trapz(par.xx(1:length(dev.krad), p_CEx(1:length(dev.krad))));

                    end
            end

            figure(3014)
            hold on
            for i = 1:length(solex_Voc.parval1)
                semilogy(solex_Voc.parval2, n_CE(i,:))
            end
            hold off
            xlabel('Light intensity [Suns]')
            ylabel('integrated electron density [cm-2]')
            hold off



        end



        function plotJV(parexsol, par1logical, par2logical)

            figure(3015)

            if length(par1logical) > length(parexsol.parval1)
                par1logical = par1logical(1:length(parexsol.parval1));
            end

            if length(par2logical) > length(parexsol.parval2)
                par2logical = par2logical(1:length(parexsol.parval2));
            end

            for i=1:length(parexsol.parval1)
                if par1logical(i) == 1
                    for j = 1:length(parexsol.parval2)
                        if par2logical(j) == 1
                            plot(squeeze(parexsol.Vapp_f(i,j,:)), -squeeze(parexsol.J_f(i,j,:)))
                            hold on
                        end
                    end
                end
            end
            xlabel('Applied Voltage [V]')
            ylabel('Current density [Acm-2]')
            ylim([-10e-3, 30e-3])
            xlim([0,1.1])
            hold off
        end

        function plotVocstable(parexsol)

            offset = parexsol.parval2-parexsol.par_base.IP(1);

            figure(3001)
            surf(offset, parexsol.parval1, parexsol.stats.Voc_stable)
            s1 = gca;
            ylabel('Mobile ion density [cm-3]')
            %ylabel('p-type SRH time constant [s]')
            xlabel('p-type VB-Fermi level offset [eV]')
            zlabel('Voc F scan [V]')
            xlim([offset(1), offset(end)]);
            ylim([parexsol.parval1(1), parexsol.parval1(end)])
            set(s1,'YScale','log');
            shading interp
            colorbar
            caxis([0.75, 0.95])
            %caxis([1.05, 1.15])
        end


        function plotJscF(parexsol)

            offset = parexsol.parval2-parexsol.par_base.IP(1);

            figure(3001)
            surf(offset, parexsol.parval1, parexsol.stats.Jsc_f)
            s1 = gca;
            %ylabel('Mobile ion density [cm-3]')
            ylabel('p-type SRH time constant [s]')
            xlabel('p-type VB-Fermi level offset [eV]')
            zlabel('Jsc F scan [Acm-2]')
            xlim([offset(1), offset(end)]);
            ylim([parexsol.parval1(1), parexsol.parval1(end)])
            set(s1,'YScale','log');
            shading interp
            colorbar
            %caxis([0.75, 0.95])
            %caxis([1.05, 1.15])
        end


    end

end
