classdef explore
    
    methods (Static)
        
        function exsol = explore2par(par_base, parnames, parvalues, JVpnts)
            % EXPLORE2PAR is used to explore 2 different parameters using the
            % parallel computing toolbox. The code is likely to require
            % modification for individual parameters owing to variable dependencies.
            % PAR_BASE is the base parameter set
            % PARNAMES is a cell array with the parameter names in - check these
            % carefully to avoid heartache later
            % PARVALUES is matrix with the parameter value ranges e.g.
            % Note: Changes in device dimensions, DCELL must be the FIRST
            % of the two parameters in PARNAMES and associated values in
            % PARVALUES
            tic
            disp('Starting parameter exploration');
            disp(['Parameter 1: ', parnames(1)]);
            disp(['Parameter 2: ', parnames(2)]);
            parval1 = cell2mat(parvalues(1));
            parval2 = cell2mat(parvalues(2));
            str1 = char(parnames(1));
            str2 = char(parnames(2));
            errorlog = zeros(length(parval1), length(parval2));

            j = 1;
            parfor i = 1:length(parval1)
                
                par = par_base;
                par.Ana = 0;
                par = explore.helper(par, str1, parval1(i));
                
                if strmatch('dcell', parnames(1)) ~= 0
                    % sets PCELL at the required position to provide a point
                    % density of one point per nanometer for changes in
                    % thickness
                    layerpoints = round(parval1(i)*1e7);
                    par = explore.helper(par, ['p', str1(2:end)], layerpoints);
                end
                
                % Refresh device
                par = pc.refresh(par);
                
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
                errortemp = zeros(1,length(parval2));
                % Stabilised profiles - Not currently working
%                 n_f = zeros(length(parval1), length(parval2), 10000);
%                 p_f = zeros(length(parval1), length(parval2), 10000);
%                 a_f = zeros(length(parval1), length(parval2), 10000);
%                 V_f = zeros(length(parval1), length(parval2), 10000);
%                 x = zeros(length(parval1), length(parval2), 10000);
            
                Voc_stable = zeros(length(parval2), JVpnts);
                PLint = zeros(length(parval2), JVpnts);
                Vapp_f = zeros(1, JVpnts);
                J_f = zeros(length(parval2), JVpnts);
                Vapp_r = zeros(1, JVpnts);
                J_r = zeros(length(parval2), JVpnts);
                
                for j = 1:length(parval2)
                    try
                        disp(['Run no. ', num2str((i-1)*length(parval2) + j), ', ', str1 ,' = ', num2str(parval1(i)), ' , ',str2,' = ', num2str(parval2(j))]);
                        
                        % If the second parameter name is intensity then set
                        % INT using PARVAL2 else set the appropriate
                        % parameter to the required value in par and Int is
                        % assumed to be 1
                        if strmatch('Int', parnames(2)) == 1
                            Int = parval2(j);
                        else
                            par = explore.helper(par, str2, parval2(j));
                            Int = 1;
                            % Refresh device
                            par = pc.refresh(par);
                        end
                        
                        % Obtain equilibrium solution - could be calculated
                        % in the i loop to avoid recalculation but depends
                        % on the variables being altered.
                        soleq = equilibrate(par);
                        
                        % Perform JV from 0-1.3 V at 50 mVs-1
                        JV = doJV(soleq.ion, 50e-3, JVpnts, Int, 1, 0, 1.3, 2);
                        
                        % Get J-V stats
                        stats = dfana.JVstats(JV);
                        
                        % Write stats into temporary variables
                        Vapp_f(j,:) = dfana.calcVapp(JV.ill.f);
                        J_f(j,:) = explore.getJtot(JV.ill.f);
                        Vapp_r(j,:) = dfana.calcVapp(JV.ill.r);
                        J_r(j,:) = explore.getJtot(JV.ill.r);
                        Voc_f(j) = stats.Voc_f;
                        Voc_r(j) = stats.Voc_r;
                        Jsc_f(j) = stats.Jsc_f;
                        Jsc_r(j) = stats.Jsc_r;
                        mpp_f(j) = stats.mpp_f;
                        mpp_r(j) = stats.mpp_r;
                        FF_f(j) = stats.FF_f;
                        FF_r(j) = stats.FF_r;
                        
                        % For steady-state solution
                        % lighton_Rs(sol_ini, int1, stable_time, mobseti, Rs, pnts)
                        sol_ill = lighton_Rs(soleq.ion, Int, -10, 1, 1e6, JVpnts);
                        
                        % Write steady-state solutions into temporary
                        % variables
                        Voc_stable(j,:) = dfana.Voct(sol_ill);
                        PLint(j,:) = dfana.PLt(sol_ill);
                        n_av(j) = mean(sol_ill.u(end, par.pcum(2):par.pcum(5),1));
                        p_av(j) = mean(sol_ill.u(end, par.pcum(2):par.pcum(5),2));
%                         n_f = explore.writevar(n_f, i, j, par.xx, sol_ill.u(end,:,1));
%                         p_f = explore.writevar(p_f, i, j, par.xx, sol_ill.u(end,:,2));
%                         a_f = explore.writevar(a_f, i, j, par.xx, sol_ill.u(end,:,3));
%                         V_f = explore.writevar(V_f, i, j, par.xx, sol_ill.u(end,:,4));
%                         x = explore.writevar(x, i, j, par.xx, sol_ill.x);
                        errortemp(j) = 0;
                    catch
                        warning(['DRIFTFUSION FAILURE: Run no. ', num2str((i-1)*length(parval2) + j), ', ', str1, '= ',num2str(parval1(i)), ', ', str2, '= ', num2str(parval2(j))]);
                        errortemp(j) = 1;
                    end
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
%                 GG = n_f;
%                 HH = p_f;
%                 II = a_f;
%                 JJ = V_f;
%                 KK = x;
                errorlog(i,:) = errortemp;
                
            end
            
            % Store solutions in output struct
            exsol.stats.Voc_f = A;
            exsol.stats.Voc_r = B;
            exsol.stats.Jsc_f = C;
            exsol.stats.Jsc_r = D;
            exsol.stats.mpp_f = E;
            exsol.stats.mpp_r = F;
            exsol.stats.FF_f = G;
            exsol.stats.FF_r = H;
            exsol.Vapp_f = AA;
            exsol.J_f = BB;
            exsol.Vapp_r = CC;
            exsol.J_r = DD;
            exsol.stats.Voc_stable = J;
            exsol.stats.PLint = K;
            exsol.parnames = parnames;
            exsol.parvalues = parvalues;
            exsol.parval1 = parval1;
            exsol.parval2 = parval2;
            exsol.par_base = par_base;
            exsol.stats.n_av = EE;
            exsol.stats.p_av = FF;
            exsol.n_f = GG;
            exsol.p_f = HH;
            exsol.a_f = II;
            exsol.V_f = JJ;
            exsol.x = KK;
            exsol.errorlog = errorlog;
            
            toc
            
        end
        
        function Jtot = getJtot(sol)
           % Extract the total current for current structure
           [~,J] = dfana.calcJ(sol);
           Jtot = J.tot(:,end);
        end
        
        function var = writevar(var, i, j, xx, arr)
            % EXPLORE.WRITEVAR writes array ARR to variable VAR using variable length
            % assignment, which is not usually allowed in PARFOR loops. This
            % allows the solution with different length vectors to be stored.
            var(i, j, 1:length(xx)) = arr;
        end
        
        function par = helper(par, parname, parvalue)
            % takes parameter set and sets parname to parvalue- workaround for parallel
            % computing loops
            eval(['par.',parname,'= parvalue;']);
        end
        
        function plotPL(exsol)
            figure(100)
            s1 = surf(exsol.parval1, exsol.parval2, exsol.stats.PLint);
            ylabel(exsol.parnames(1))
            xlabel(exsol.parnames(2))
            set(s1,'YScale','log');
            zlabel('PL intensity [cm-2s-1]')
            shading interp
            colorbar
        end
        
        function plotsurf(exsol, yproperty, xlogon, ylogon, zlogon)
            % Plots a surface with a
            try
                eval(['y = exsol.stats.', yproperty, ';'])
            catch
                error('YPROPERTY is not a property contained in exsol.STATS- please check YPROPERTY string')
            end
            
            figure(101)
            surf(exsol.parval2, exsol.parval1, exsol.stats.Voc_f)
            s1 = gca;
            ylabel(exsol.parnames{1,1})
            xlabel(exsol.parnames{1,2})
            zlabel(yproperty)
            xlim([exsol.parval2(1), exsol.parval2(end)])
            ylim([exsol.parval1(1), exsol.parval1(end)])
            
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
            shading interp
            colorbar
            cb = colorbar();
            if zlogon
                cb.Ruler.Scale = 'log';
                cb.Ruler.MinorTick = 'on';
            end
        end
        
        function plotstat_2D_parval1(exsol, yproperty, logx, logy)
            % YPROPERTY is string with the property name. Properties must
            % be one of those contained in the exsol.STATS structure
            eval(['y = exsol.stats.', yproperty, ';'])
            for i=1:length(exsol.parval1)
                figure(103)
                if logx == 0 && logy == 0
                    plot(exsol.parval2, y(i,:))
                elseif logx == 1 && logy == 0
                    semilogx(exsol.parval2, y(i,:))
                elseif logx == 0 && logy == 1
                    semilogy(exsol.parval2, y(i,:))
                elseif logx == 1 && logy == 1
                    loglog(exsol.parval2, y(i,:))
                end
                hold on
            end
            xlabel(exsol.parnames{1,2})
            ylabel(yproperty)
            legstr = (num2str((exsol.parval1'+60e-7)*1e7));
            legend(legstr)
            xlim([min(exsol.parval2), max(exsol.parval2)])
            hold off
        end
        
        function plotstat_2D_parval2(exsol, yproperty, logx, logy)
            % YPROPERTY is string with the property name. Properties must
            % be one of those contained in the exsol.STATS structure
            figure(104)
            
            eval(['y = exsol.stats.', yproperty,';'])
            
            if strmatch('Voc_stable', yproperty) ==1
                y = y(:,:,end);
            end
            
            if strmatch('Jsc_r', yproperty) ==1
                y = y;
            end
            
            x = (exsol.parval1'+ 60e-7)*1e7;
            
            for i=1:length(exsol.parval2)
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
            %xlabel(exsol.parnames{1,1})
            xlabel('Active layer thickness [nm]')
            ylabel(yproperty)
            legstr = (num2str(exsol.parval2'));
            legend(legstr)
            xlim([min(x), max(x)])
            
            hold off
        end
        
        function plotfinalELx(exsol)
            % plots final energy level diagrams
            figure(105)
            for i=1:length(parval1)
                for j = 1:length(parval2)
                    % Rebuild solutions
                    sol.u(1,:,1) = exsol.nf(i, j, :);
                    sol.u(1,:,2) = exsol.pf(i, j, :);
                    sol.u(1,:,3) = exsol.af(i, j, :);
                    sol.u(1,:,4) = exsol.Vf(i, j, :);
                    sol.t = 0;
                    sol.x = exsol.x(i,j,:);
                    sol.par = exsol.par_base;
                    
                end
            end
        end
        
        function plotprof_2D(exsol, yproperty, par1logical, par2logical, logx,logy)
            % PLOTPROF_2D plots one dimensional profiles of YPROPERTY using
            % PAR1LOGICAL and PAR2LOGICAL to determine which solutions to
            % plot.
            % In this case YPROPERTY must be a three dimensional vector with
            % dimensions [length(parval1), length(parval2), length(x)]
            % contained within EXSOL
            par = exsol.par_base;
            
            eval(['y = exsol.', yproperty,';']);
            if length(par1logical) > length(exsol.parval1)
                par1logical = par1logical(1:length(exsol.parval1));
            end
            
            if length(par2logical) > length(exsol.parval2)
                par2logical = par2logical(1:length(exsol.parval2));
            end
            
            % Replace zeros with NaNs
            y(y==0) = NaN;
            exsol.x(exsol.x==0) = NaN;
            
            if strmatch(yproperty,'a_f')
                y = y-exsol.par_base.Nani(1);
            end
            
            parval1 = cell2mat(exsol.parvalues(1));
            parval2 = cell2mat(exsol.parvalues(2));
            str1 = char(exsol.parnames(1));
            str2 = char(exsol.parnames(2));
            
            figure(106)
            for i=1:length(exsol.parval1)
                if par1logical(i) == 1
                    
                    if strmatch('dcell', str1) ~= 0
                        % sets PCELL at the required position to provide a point
                        % density of one point per nanometer for changes in
                        % thickness
                        layerpoints = round(parval1(i)*1e7);
                        par = explore.helper(par, ['p', str1(2:end)], layerpoints);
                    end
                    
                    % Refresh device
                    par = pc.refresh(par);
                    
                    for j = 1:length(exsol.parval2)
                        if par2logical(j) == 1
                            % Rebuild solutions
                            if logx
                                semilogx(squeeze(exsol.x(i, j, :).*1e7), squeeze(y(i, j, :)));
                            elseif logy
                                semilogy(squeeze(exsol.x(i, j, :).*1e7), squeeze(y(i, j, :)));
                            elseif logx ==1 && logy ==1
                                loglog(squeeze(exsol.x(i, j, :).*1e7), squeeze(y(i, j, :)));
                            else
                                plot(squeeze(exsol.x(i, j, :).*1e7), squeeze(y(i, j, :)));
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
        
        function plotU(exsol, par1logical, par2logical,logx,logy)
            
            if length(par1logical) > length(exsol.parval1)
                par1logical = par1logical(1:length(exsol.parval1));
            end
            
            if length(par2logical) > length(exsol.parval2)
                par2logical = par2logical(1:length(exsol.parval2));
            end
            
            UHTL = zeros(length(par1logical),length(par2logical));
            UHTLint = zeros(length(par1logical),length(par2logical));
            Ubulk = zeros(length(par1logical),length(par2logical));
            UETLint = zeros(length(par1logical),length(par2logical));
            UETL = zeros(length(par1logical),length(par2logical));
            Utotal = zeros(length(par1logical),length(par2logical));
            
            figure(107)
            par = exsol.par_base;
            
            for i=1:length(exsol.parval1)
                if par1logical(i) == 1
                    for j = 1:length(exsol.parval2)
                        if par2logical(j) == 1
                            % Rebuild solutions
                            par = explore.helper(par, exsol.parnames{1,1}, exsol.parval1(i));
                            
                            if strmatch('dcell(1,4)', exsol.parnames(1)) ~= 0
                                pcontact = round(exsol.parval1(i)*1e7);
                                par.pcell(1,4) = pcontact*1;
                            end
                            
                            % Refresh device
                            par = pc.refresh(par);
                            
                            dev = par.dev;
                            
                            n = squeeze(exsol.n_f(i,j,:));
                            p = squeeze(exsol.p_f(i,j,:));
                            n=n';
                            p=p';
                            n = n(1:length(dev.krad));
                            p = p(1:length(dev.krad));
                            
                            % Recombination
                            Ubtb = dev.krad.*(n.*p - dev.ni.^2);
                            
                            Usrh = ((n.*p - dev.ni.^2)./((dev.taun.*(p+dev.pt)) + (dev.taup.*(n+dev.nt))));
                            
                            U = Ubtb + Usrh;
                            
                            xplot = squeeze(exsol.x(i, j, :).*1e7);
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
            
            figure(108)
            xlabel('Position [nm]')
            ylabel('U [cm-3s-1]')
            hold off
            
            for j=1:length(exsol.parval2)
                if par2logical(j) == 1
                    figure(114)
                    plot(exsol.parval1*1e7, UHTL(:,j),...
                        exsol.parval1*1e7, UHTLint(:,j),...
                        exsol.parval1*1e7, Ubulk(:,j),...
                        exsol.parval1*1e7, UETLint(:,j),...
                        exsol.parval1*1e7, UETL(:,j),...
                        exsol.parval1*1e7, Utotal(:,j));
                    %hold on
                end
            end
            xlabel('Active layer thickness [nm]')
            ylabel('Recombination rate [cm-3s-1]')
            legend('HTL', 'HTL interface', 'Bulk perovskite', 'ETL interface', 'ETL','Total device')
            hold off
            
        end
        
        function plotCE(exsol_Voc, exsol_eq, xlogon, ylogon, zlogon, normalise)
            % plotting function for plotting carrier densities as a
            % function of Voc
            
            % exsol_Voc = EXPLORE open circuit solutions
            % exsol_eq = EXPLORE equilibrium solutions
            % The carrier densities at equilibrium are subtracted from
            % those at Voc to obtain the 'extracted' electronic carrier
            % densities
            % This version requires that the equilibrium solutions are
            % found using a separate run of EXPLORE.EXPLORE2PAR with only a
            % single J element
            
            par = exsol_Voc.par_base;
            parval1 = cell2mat(exsol_Voc.parvalues(1));
            parval2 = cell2mat(exsol_Voc.parvalues(2));
            str1 = char(exsol_Voc.parnames(1));
            str2 = char(exsol_Voc.parnames(2));
            
            n_CE = zeros(length(exsol_Voc.parval1), length(exsol_Voc.parval2));
            p_CE = zeros(length(exsol_Voc.parval1), length(exsol_Voc.parval2));
            
            for i=1:length(exsol_Voc.parval1)   % d_active
                % Rebuild solutions
                par = explore.helper(par, exsol_Voc.parnames{1,1}, exsol_Voc.parval1(i));
                
                if strmatch('dcell', str1) ~= 0
                    % sets PCELL at the required position to provide a point
                    % density of one point per nanometer for changes in
                    % thickness
                    layerpoints = round(parval1(i)*1e7);
                    par = explore.helper(par, ['p', str1(2:end)], layerpoints);
                end
                
                % Refresh device
                par = pc.refresh(par);
                x = par.xx;
                dev = par.dev;
                
                n_CEx = zeros(1,length(par.xx));
                p_CEx = zeros(1,length(par.xx));
                
                for j = 1:length(exsol_Voc.parval2) % Light intensity
                        
                    n_Voc = squeeze(exsol_Voc.n_f(i,j,1:length(par.xx)))';
                    p_Voc = squeeze(exsol_Voc.p_f(i,j,1:length(par.xx)))';
                    n_eq = squeeze(exsol_eq.n_f(i,1:length(par.xx)));
                    p_eq = squeeze(exsol_eq.p_f(i,1:length(par.xx)));
                    
                    n_CEx = n_Voc - n_eq;
                    p_CEx = p_Voc - p_eq;
                    
                    n_CE(i,j) = trapz(par.xx, n_CEx);
                    p_CE(i,j) = trapz(par.xx, p_CEx);
                    
                    if normali
                end
            end
            
            figure(109)
            hold on
            for i = 1:length(exsol_Voc.parval1)
                loglog(exsol_Voc.parval2, n_CE(i,:))
            end
            hold off
            xlabel('Light intensity [Suns]')
            ylabel('integrated electron density [cm-2]')
            hold off
            
            figure(110)
            d_active = round(exsol_Voc.parval1*1e7)+60;
            surf(exsol_Voc.parval2, d_active , n_CE)
            s1 = gca;
%             ylabel(exsol_Voc.parnames{1,1})
%             xlabel(exsol_Voc.parnames{1,2})
            xlabel('Light intensity [Sun]')
            ylabel('Active layer thickness [nm]')
            zlabel('integrated electron density [cm-2]')
            xlim([exsol_Voc.parval2(1), exsol_Voc.parval2(end)])
            ylim([d_active(1),d_active(end)])
            
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
            shading interp
            colorbar
            cb = colorbar();
            if zlogon
                cb.Ruler.Scale = 'log';
                cb.Ruler.MinorTick = 'on';
            end
            
            figure(111)
            d_active = round(exsol_Voc.parval1*1e7)+60;
            surf(exsol_Voc.parval2, d_active, p_CE)
            s1 = gca;
%             ylabel(exsol_Voc.parnames{1,1})
%             xlabel(exsol_Voc.parnames{1,2})
            xlabel('Light intensity [Sun]')
            ylabel('Active layer thickness [nm]')
            zlabel('integrated hole density [cm-2]')
            xlim([exsol_Voc.parval2(1), exsol_Voc.parval2(end)])
            ylim([d_active(1),d_active(end)])
            
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
            shading interp
            colorbar
            cb = colorbar();
            if zlogon
                cb.Ruler.Scale = 'log';
                cb.Ruler.MinorTick = 'on';
            end
            end
        end
        
        function plotJV(exsol, par1logical, par2logical)
            figure(112)
            
            if length(par1logical) > length(exsol.parval1)
                par1logical = par1logical(1:length(exsol.parval1));
            end
            
            if length(par2logical) > length(exsol.parval2)
                par2logical = par2logical(1:length(exsol.parval2));
            end
            
            for i=1:length(exsol.parval1)
                if par1logical(i) == 1
                    for j = 1:length(exsol.parval2)
                        if par2logical(j) == 1
                            plot(squeeze(exsol.Vapp_f(i,j,:)), squeeze(exsol.J_f(i,j,:)))
                            hold on
                        end
                    end
                end
            end
            xlabel('Applied Voltage [V]')
            ylabel('Current density [Acm-2]')
            ylim([-30e-3, 10e-3])
            hold off
        end
        
    end

end
