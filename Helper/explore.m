classdef explore
    
    methods (Static)
        
        function parexsol = explore2par(par_base, parnames, parvalues, JVpnts)
            % EXPLOREPAR is used to explore 2 different parameters using a parallel pool.
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
                
                par = par_base;
                par.Ana = 0;
                par = explore.helper(par, str1, parval1(i));
                
                if strmatch('dcell(1,4)', parnames(1)) ~= 0
                    pcontact = round(parval1(i)*1e7);
                    par.pcell(1,4) = pcontact*1;
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
                    disp(['Run no. ', num2str(runN), ', d_active = ', num2str(1e7*(parval1(i)+60e-7)), ' nm, Intensity = ', num2str(parval2(j))]);
                    
                    %par = explore.helper(par, str2, parval2(j));
                    
                    % Variable 2 will be light intensity
                    JV = doJV(soleq.ion, 50e-3, JVpnts, parval2(j), 1, 0, 1.3, 2);
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
                    
                    sol_ill = lighton(soleq.ion, parval2(j), 10, 1, 1e6, JVpnts)
                    
                    % For PL
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
                J(:,:,i) = Voc_stable;
                K(:,:,i) = PLint;
                AA(:,:,i) = Vapp_f;
                BB(:,:,i) = J_f;
                CC(:,:,i) = Vapp_r;
                DD(:,:,i) = J_r;
                EE(i,:) = n_av;
                FF(i,:) = p_av;
                GG(:,:,i) = n_f;
                HH(:,:,i) = p_f;
                II(:,:,i) = a_f;
                JJ(:,:,i) = V_f;
                KK(:,:,i) = x;
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
            parexsol.Vapp_f = CC;
            parexsol.J_f = DD;
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
            
            save('ws_explore.mat')  % Save workspace
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
            eval(['y = parexsol.stats.', yproperty])
            x = (parexsol.parval1'+ 60e-7)*1e7;
            for i=1:length(parexsol.parval2)
                figure(3010)
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
