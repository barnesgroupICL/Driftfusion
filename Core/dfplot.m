classdef dfplot
    % Plotting class - contains methods for plotting
    % DFPLOT.ELX = Energy level diagrams and charge carrier densities
    % DFPLOT.JT = Currents as a function of time
    % DFPLOT.JX = Total currents as a function of position
    % DFPLOT.JDDX = Drift and diffusion currents as a function of position
    % Plotting functions that are a funciton of position can accept a time
    % array as the second argument- the procedure will loop and plot the
    % solution at multiple times.
    % The third optional argument defines the x-range.
    
    methods (Static)
        
        function ELxnpxacx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            
            figure(1);
            subplot(3,1,1);
            dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'\it{E}_{fn}', '\it{E}_{fp}', '\it{E}_{CB}', '\it{E}_{VB}'}, {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)
            
            subplot(3,1,2);
            dfplot.x2d(sol, x, {n, p}, {'\it{n}', '\it{p}'}, {'-', '-.'}, 'El. carrier density [cm-3]', tarr, xrange, 0, 1)
            
            subplot(3,1,3);
            dfplot.x2d(sol, x, {c, a}, {'\it{c}', '\it{a}'}, {'-', '-.'}, 'Ionic carrier density [cm-3]', tarr, xrange, 0, 0)
        end

        function ELxnpx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            
            figure(101);
            subplot(2,1,1);
            dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'\it{E}_{fn}', '\it{E}_{fp}', '\it{E}_{CB}', '\it{E}_{VB}'}, {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)
            
            subplot(2,1,2);
            dfplot.x2d(sol, x, {n, p}, {'\it{n}', '\it{p}'}, {'-', '-.'}, 'El. carrier density [cm-3]', tarr, xrange, 0, 1)
            
        end
        
        function acnpFx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            [FV, Frho] = dfana.calcF(sol);
            NAmat = repmat(dev.NA, length(t), 1);
            NDmat = repmat(dev.ND, length(t), 1);
            
            figure(26);
            subplot(3,1,1);
            dfplot.x2d(sol, x, {a, c}, {'a', 'c'}, {'-.', '-'}, 'Ionic carrier density [cm-3]', tarr, xrange, 0, 0)
            
            subplot(3,1,2);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-.'}, 'Electronic carrier density [cm-3]', tarr, xrange, 0, 1)
            
            subplot(3,1,3);
            dfplot.x2d(sol, x, {FV}, {''}, {'-'}, 'Electric field [Vcm-1]', tarr, xrange, 0, 0)
        end

        function acNANDnpFx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            [FV, Frho] = dfana.calcF(sol);
            NAmat = repmat(dev.NA, length(t), 1);
            NDmat = repmat(dev.ND, length(t), 1);
            
            figure(26);
            subplot(3,1,1);
            dfplot.x2d(sol, x, {a + NAmat, c + NDmat}, {'NA', 'ND'}, {'-', '-.'}, 'Dopant density [cm-3]', tarr, xrange, 0, 0)
            
            subplot(3,1,2);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-.'}, 'Electronic carrier density [cm-3]', tarr, xrange, 0, 1)
            
            subplot(3,1,3);
            dfplot.x2d(sol, x, {FV}, {''}, {'-'}, 'Electric field [Vcm-1]', tarr, xrange, 0, 0)
        end
        
        function acnpx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            NAmat = repmat(dev.NA, length(t), 1);
            NDmat = repmat(dev.ND, length(t), 1);
            
            figure(26);
            subplot(2,1,1);
            dfplot.x2d(sol, x, {a, c}, {'a', 'c'}, {'-', '-', '--', '--'}, 'Ionic carrier density [cm-3]', tarr, xrange, 0, 0)
            
            subplot(2,1,2);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-'}, 'Electronic carrier density [cm-3]', tarr, xrange, 0, 1)
        end

        function rhoacnpx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            NAmat = repmat(dev.NA, length(t), 1);
            NDmat = repmat(dev.ND, length(t), 1);
            
            figure(26);
            subplot(2,1,1);
            dfplot.x2d(sol, x, {c-a+NDmat-NAmat}, {''}, {'-'}, 'Ionic space charge density [cm-3]', tarr, xrange, 0, 0)
            
            subplot(2,1,2);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-'}, 'Electronic carrier density [cm-3]', tarr, xrange, 0, 1)
        end
        
        function rhoacnpFx2(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            [FV, Frho] = dfana.calcF(sol);
            NAmat = repmat(dev.NA, length(t), 1);
            NDmat = repmat(dev.ND, length(t), 1);
            
            figure(261);
            subplot(3,1,1);
            dfplot.x2d(sol, x, {c-a, NDmat-NAmat}, {'Mobile ions', 'Static dopants'}, {'-', '-.'}, 'Ionic space charge density [cm-3]', tarr, xrange, 0, 0)
            
            subplot(3,1,2);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-.'}, 'Electronic carrier density [cm-3]', tarr, xrange, 0, 1)
 
            subplot(3,1,3);
            dfplot.x2d(sol, x, {FV}, {''}, {'-'}, 'Electric field [Vcm-1]', tarr, xrange, 0, 0)
        end
        
        function rhoacnpVx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            [FV, Frho] = dfana.calcF(sol);
            NAmat = repmat(dev.NA, length(t), 1);
            NDmat = repmat(dev.ND, length(t), 1);
            
            figure(261);
            subplot(3,1,1);
            dfplot.x2d(sol, x, {c-a, NDmat-NAmat}, {'Mobile ions', 'Static dopants'}, {'-', '-.'}, 'Ionic space charge density [cm-3]', tarr, xrange, 0, 0)
            
            subplot(3,1,2);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-.'}, 'Electronic carrier density [cm-3]', tarr, xrange, 0, 1)
 
            subplot(3,1,3);
            dfplot.x2d(sol, x, {V}, {''}, {'-'}, 'Electrostatic potential [V]', tarr, xrange, 0, 0)
        end
        
        function Jt(sol, xpos)
            % Currents as a function of time
            % POS = the readout position
            t = sol.t;
            [J, j, xmesh] = dfana.calcJ(sol);
            ppos = getpointpos(xpos, xmesh);
            
            figure(2);
            plot(t, J.n(:, ppos),t, J.p(:, ppos),t, J.a(:, ppos),t, J.c(:, ppos), t, J.disp(:,ppos), t, J.tot(:, ppos));
            legend('Jn', 'Jp', 'Ja', 'Jc', 'Jdisp', 'Jtotal')
            xlabel('Time [s]');
            ylabel('J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end

        function Jtott(sol, xpos)
            % Currents as a function of time
            % POS = the readout position
            t = sol.t;
            [J, j, xmesh] = dfana.calcJ(sol);
            ppos = getpointpos(xpos, xmesh);
            
            figure(2);
            plot(t, J.tot(:, ppos));
            %legend('Jn', 'Jp', 'Ja', 'Jc', 'Jdisp', 'Jtotal')
            xlabel('Time [s]');
            ylabel('J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end
        
        function Jx(varargin)
            % Plots the currents
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [J, j, x] = dfana.calcJ(sol);
            
            figure(3);
            dfplot.x2d(sol, x, {J.n, J.p, J.a, J.c, J.disp, J.tot},...
                {'Jn', 'Jp', 'Ja', 'Jc', 'Jdisp', 'Jtot'}, {'-','-','-','-','-','-'},...
                'Current density [Acm-2]', tarr, xrange, 0, 0);
        end
        
        function jx(varargin)
            % Plots the carrier fluxes
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [J, j, x] = dfana.calcJ(sol);
            
            figure(301);
            dfplot.x2d(sol, par.x_ihalf, {j.n, j.p, j.a, j.c, j.disp},{'jn', 'jp', 'ja', 'jc', 'jdisp'},...
                {'-','-','-','-','-'}, 'Current density [Acm-2]', tarr, xrange, 0, 0);
        end
        
        function JV(JV, option)
            % JV - a solution from doJV
            % OPTION - 1 = dark only, 2 = light only, 3 = dark & light
            % JV is a structure containing dark and illuminated JVs
            
            if option == 1 || option == 3
                J.dk.f = dfana.calcJ(JV.dk.f);
                Vapp.dk.f = dfana.calcVapp(JV.dk.f);
                J.dk.r = dfana.calcJ(JV.dk.r);
                Vapp.dk.r = dfana.calcVapp(JV.dk.r);
                
                figure(4)
                plot(-Vapp.dk.f, J.dk.f.tot(:,end), '--', -Vapp.dk.r, J.dk.r.tot(:,end));
                hold on
            end
            
            if option == 2 || option == 3
                
                J.ill.f = dfana.calcJ(JV.ill.f);
                Vapp.ill.f = dfana.calcVapp(JV.ill.f);
                J.ill.r = dfana.calcJ(JV.ill.r);
                Vapp.ill.r = dfana.calcVapp(JV.ill.r);
                
                figure(4)
                plot(-Vapp.ill.f, J.ill.f.tot(:,end),'--')%, 'Color', [0, 0.4470, 0.7410]);
                hold on
                plot(-Vapp.ill.r, J.ill.r.tot(:,end));%,'Color', [0, 0.4470, 0.7410]);
            end
            
            figure(4)
            %ylim([-30e-3, 10e-3]);
            xlabel('Applied voltage [V]')
            ylabel('Current density [Acm-2]');
            hold off
            
        end
        
        function Jddx(varargin)
            % figure(5)
            % drift and diffusion currents as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [~, Jdd, x] = dfana.Jddxt(sol);
            
            figure(301);
            dfplot.x2d(sol, x, {Jdd.ndiff, Jdd.ndrift, Jdd.pdiff, Jdd.pdrift,...
                Jdd.adiff, Jdd.adrift, Jdd.cdiff, Jdd.cdrift},...
                {'Jn,diff', 'Jn,drift', 'Jp,diff', 'Jp,drift', 'Ja,diff', 'Ja,drift', 'Jc,diff', 'Jc,drift'},...
                {'.','-','.','-','.','-','.','-'},'Current density [Acm-2]', tarr, xrange, 0, 0);
        end
        
        function Voct(sol)
            Voc = dfana.calcVQFL(sol);
            figure(6)
            plot(sol.t, Voc)
            xlabel('Time [s]')
            ylabel('Voc [V]')
        end
        
        function PLt(sol)
            PL = dfana.PLt(sol);
            figure(7)
            plot(sol.t, PL)
            xlabel('Time [s]')
            ylabel('PL [cm-2s-1]')
        end
        
        function Vappt(sol)
            % Difference in potential between the left and right boundary
            figure(8)
            plot(sol.t, -(sol.u(:,end,4)-sol.u(:,1,4)-sol.par.Vbi));
            xlabel('Time [s]')
            ylabel('Vapp [V]')
        end
        
        function JVapp(sol, xpos)
            % Obtain point position from x position
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            
            J = dfana.calcJ(sol);
            Vapp = -(sol.u(:,end,4)-sol.u(:,1,4)-sol.par.Vbi);
            
            figure(9)
            plot(Vapp, J.n(:, ppos),Vapp, J.p(:, ppos),Vapp, J.a(:, ppos),Vapp, J.disp(:,ppos), Vapp, J.tot(:, ppos));
            legend('Jn', 'Jp', 'Ja', 'Jdisp', 'Jtotal')
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end
        
        function JtotVapp(sol, xpos)
            % Obtain point position from x position
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            
            J = dfana.calcJ(sol);
            Vapp = -(sol.u(:,end,4)-sol.u(:,1,4)-sol.par.Vbi);
            
            figure(91)
            plot(Vapp, J.tot(:, ppos));
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end
        
        
        function logJVapp(sol, xpos)
            % plot the log of the mod J
            
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            
            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol);
            
            figure(10)
            semilogy(Vapp, abs(J.tot(:,ppos)), Vapp, abs(J.n(:,ppos)), Vapp, abs(J.p(:,ppos)), Vapp, abs(J.a(:,ppos)),Vapp, abs(J.c(:,ppos)), Vapp, abs(J.disp(:,ppos)));
            xlabel('Vapp [V]');
            ylabel('|J| [A cm^{-2}]');
            legend('Jtot', 'Jn', 'Jp', 'Ja', 'Jc', 'Jdisp')
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end
        
        function logJVapp3D(sol, xpos, ylogon)
            
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            
            t = sol.t;
            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol)';
            Jtot=J.tot(:, ppos);
            
            figure(11)
            surface('XData', [Vapp Vapp],             ... % N.B.  XYZC Data must have at least 2 cols
                'YData', [abs(Jtot) abs(Jtot)],             ...
                'ZData', [t' t'], ...
                'CData', [t' t'],             ...
                'FaceColor', 'none',        ...
                'EdgeColor', 'interp',      ...
                'Marker', 'none','LineWidth',1);
            s1 = gca;
            xlabel('Vapp [V]');
            ylabel('|J| [A cm^{-2}]');
            
            if ylogon
                set(s1,'YScale','log');
            else
                set(s1,'YScale','linear');
            end
            hold off
        end
        
        function xmesh(sol)
            figure(11)
            plot(sol.x)
            xlabel('Position [cm]')
            ylabel('Point')
        end
        
        function Vx(varargin)
            % Electrostatic potential as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            figure(12);
            dfplot.x2d(sol, x, {V},{'V'},{'-'},'Electrostatic potential [V]', tarr, xrange, 0, 0);
        end
        
        function Fx(varargin)
            % Electrostatic potential as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            F = dfana.calcF(sol);
            
            figure(12);
            dfplot.x2d(sol, x, {F},{'F'},{'-'},'Electric Field [Vcm-1]', tarr, xrange, 0, 0);
        end
        function npx(varargin)
            % Carrier densities as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            figure(13);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-','-'},...
                'Carrier density [cm-3]', tarr, xrange, 0, 1)
        end
        
        function acx(varargin)
            % Ionic carrier densities as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            figure(14)
            dfplot.x2d(sol, x, {a,c},{'anion','cation'}, {'-','-'},...
                'Ionic carrier density [cm-3]', tarr, xrange, 0, 0);
        end
        
        function rhoacxNAxNDx(varargin)
            % Ionic carrier densities as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            NAmat = repmat(dev.NA, length(t), 1);
            NDmat = repmat(dev.ND, length(t), 1);
            
            % charge density
            rho_ion = - a + c;
            
            figure(14)
            dfplot.x2d(sol, x, {rho_ion, NAmat, NDmat},{'Cations', 'NA', 'ND'}, {'-','--','--'},...
                'Charge Density [cm-3]', tarr, xrange, 0, 0);
        end      
        
        function gx(varargin)
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [g1, g2, g] = dfana.calcg(sol);
            
            figure(15)
            dfplot.x2d(sol, par.x_ihalf, {g1, g2, g}, {'g1', 'g2', 'g total'},...
                {'-','-','-'}, 'Generation rate [cm-3s-1]', tarr, xrange, 0, 0);
        end
        
        function gxt(sol)
            % Carrier densities as a function of position
            par = sol.par;
            t = sol.t;
            [~, ~, g] = dfana.calcg(sol);
            xnm = par.x_ihalf*1e7;
            
            figure(16)
            surf(xnm, sol.t, g)
            xlabel('Position [cm]')
            ylabel('Time [s]')
            zlabel('Generation rate [cm^{-3}s^{-1}]')
        end
        
        function Ux(varargin)
            % Recombination rates as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            U = dfana.calcU(sol);
            
            figure(17)
            dfplot.x2d(sol, x, {U.btb, U.srh, U.tot},{'Ubtb', 'Usrh', 'Utot'},...
                {'-','-','-'}, 'Recombination rate [cm-3s-1]', tarr, xrange, 0, 0);
        end
        
        function JVrec(JV, option)
            % Plots recombination currents for JV
            
            % JV - a solution from doJV
            % OPTION - 1 = dark only, 2 = light only, 3 = dark & light
            % JV is a structure containing dark and illuminated JVs
            
            if option == 1 || option == 3
                J.dk.f = dfana.calcJ(JV.dk.f);
                Vapp.dk.f = dfana.calcVapp(JV.dk.f);
                J.dk.r = dfana.calcJ(JV.dk.r);
                Vapp.dk.r = dfana.calcVapp(JV.dk.r);
                
                figure(13)
                plot(Vapp.dk.f, J.dk.f.tot(:,end), '--', Vapp.dk.r, J.dk.r.tot(:,end));
                hold on
            end
            
            if option == 2 || option == 3
                solf = JV.ill.f;
                solr = JV.ill.r;
                par = solf.par;
                pcum0 = par.pcum0;
                
                J.ill.f = dfana.calcJ(JV.ill.f);
                Vapp.ill.f = dfana.calcVapp(JV.ill.f);
                J.ill.r = dfana.calcJ(JV.ill.r);
                Vapp.ill.r = dfana.calcVapp(JV.ill.r);
                
                U_f = dfana.calcU(JV.ill.f);
                Jrec_btb_f = JV.ill.f.par.e*trapz(JV.ill.f.x, U_f.btb, 2);
                Jrec_srhint_f = JV.ill.f.par.e*trapz(JV.ill.f.x(pcum0(2)+1:pcum0(3)), U_f.srh(:,pcum0(2)+1:pcum0(3)), 2)...
                    +JV.ill.f.par.e*trapz(JV.ill.f.x(pcum0(4)+1:pcum0(5)), U_f.srh(:,pcum0(4)+1:pcum0(5)), 2);
                Jrec_srhbulk_f = JV.ill.f.par.e*trapz(JV.ill.f.x(pcum0(3)+1:pcum0(4)), U_f.srh(:,pcum0(3)+1:pcum0(4)), 2);
                Jrec_tot_f = JV.ill.f.par.e*trapz(JV.ill.f.x, U_f.tot, 2);
                
                U_r = dfana.calcU(JV.ill.r);
                Jrec_btb_r = JV.ill.f.par.e*trapz(JV.ill.r.x, U_r.btb, 2);
                Jrec_srhint_r = JV.ill.r.par.e*trapz(JV.ill.r.x(pcum0(2)+1:pcum0(3)), U_r.srh(:,pcum0(2)+1:pcum0(3)), 2)...
                    +JV.ill.r.par.e*trapz(JV.ill.r.x(pcum0(4)+1:pcum0(5)), U_r.srh(:,pcum0(4)+1:pcum0(5)), 2);
                Jrec_srhbulk_r = JV.ill.r.par.e*trapz(JV.ill.r.x(pcum0(3)+1:pcum0(4)), U_r.srh(:,pcum0(3)+1:pcum0(4)), 2);
                Jrec_tot_r = JV.ill.r.par.e*trapz(JV.ill.r.x, U_r.tot, 2);
                
                cc=lines(4);
                
                figure(13)
                plot(Vapp.ill.f, J.ill.f.tot(:,end),'--', 'Color', cc(1,:));
                hold on
                plot(Vapp.ill.r, J.ill.r.tot(:,end), 'Color', cc(1,:));
                % Recombination currents
                plot(Vapp.ill.f, Jrec_btb_f,'--','Color', cc(2,:));
                plot(Vapp.ill.r, Jrec_btb_r,'Color', cc(2,:));
                plot(Vapp.ill.f, Jrec_srhint_f,'--','Color', cc(3,:));
                plot(Vapp.ill.r, Jrec_srhint_r, 'Color', cc(3,:));
                plot(Vapp.ill.f, Jrec_srhbulk_f,'--','Color', cc(4,:));
                plot(Vapp.ill.r, Jrec_srhbulk_r, 'Color', cc(4,:));
            end
            
            figure(13)
            ylim([-30e-3, 10e-3]);
            xlabel('Applied voltage [V]')
            ylabel('Current density [Acm-2]');
            hold off
            legend('Illumated for', 'Illumated rev','Jrec,btb for','Jrec,btb rev'...
                ,'Jrec,srh-int for','Jrec,srh-int rev','Jrec,srh-bulk for','Jrec,srh-bulk rev')
        end
        
        function Ft(sol, xpos)
            % Absolute field strength F as a function of time at point
            % position XPOS
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            
            F = dfana.calcF(sol);
            
            figure(14)
            plot(sol.t, F(:,ppos))
            xlabel('Time [s]')
            ylabel(['Electric Field at pos x = ', num2str(round(xpos*1e7)), 'nm [Vcm-1]'])
        end
        
        function sigmat(sol)
            % Plot the integrated space charge density [cm-2] as a function of time
            sigma = dfana.calcsigma(sol);
            t = sol.t;
            
            figure(15)
            plot(t, sigma)
            xlabel('Time [s]')
            ylabel('sigma [C cm-2]')
        end
        
        function Qt(sol, x1, x2)
            % Plot the integrated space charge density in Coulombs [Ccm-2] as a function of time
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            p1 = find(x<=x1);
            p1 = p1(end);
            p2 = find(x<=x2);
            p2 = p2(end);
            
            sigma = dfana.calcsigma(sol);
            Q = par.e*sigma(:, x1:x2);
            
            figure(17)
            plot(t, Q)
            xlabel('Time [s]')
            ylabel('Charge [C cm-2]')
            xlim([t(1), t(end)])
        end
        
        function QVapp(sol, x1, x2)
            % Integrated charge density as a function of applied voltage
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            p1 = find(x<=x1);
            p1 = p1(end);
            p2 = find(x<=x2);
            p2 = p2(end);
            
            sigma = dfana.calcsigma(sol);
            Vapp = dfana.calcVapp(sol);
            Q = par.e*sigma(:, x1:x2);
            
            figure(17)
            plot(Vapp, Q)
            xlabel('Vapp [V]')
            ylabel('Charge [C cm-2]')
            xlim([Vapp(1), Vapp(end)])
        end
        
        function rhox(varargin)
            % Volumetric charge density (rho) as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            rho = dfana.calcrho(sol);
            
            figure(19)
            dfplot.x2d(sol, x, {rho},{'\rho'},{'-'},'Charge density [cm-3]', tarr, xrange, 0, 0);
        end
        
        function rho_electronic_x(varargin)
            % Volumetric charge density (rho) as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            rho = dfana.calcrho_electronic(sol);
            
            figure(19)
            dfplot.x2d(sol, x, {rho},{'\rho'},{'-'},'Charge density [cm-3]', tarr, xrange, 0, 0);
        end
        
        function deltarhox(varargin)
            % The change in volumetric charge density (rho) as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            rho = dfana.calcrho(sol);
            deltarho = rho - rho(1,:);
            
            figure(20)
            dfplot.x2d(sol, x, {deltarho},{'\Delta \rho'},{'-'},'Delta charge density [cm-3]', tarr, xrange, 0, 0);
        end
        
        function rhoxFxVx(varargin)
            % Three panel figure:
            % Volumetric charge density (rho), Field and potential as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            rho = dfana.calcrho(sol);
            F = dfana.calcF(sol);
            
            figure(21)
            subplot(3, 1, 1)
            dfplot.x2d(sol, x, {rho},{'\rho'},{'-'}, 'Charge density [cm-3]', tarr, xrange, 0, 0);
            
            subplot(3, 1, 2)
            dfplot.x2d(sol, x, {F},{'F'},{'-'},'Electric field [Vcm-1]', tarr, xrange, 0, 0);
            
            subplot(3, 1, 3)
            dfplot.x2d(sol, x, {V},{'V'},{'-'},'Electrostatic potential [V]', tarr, xrange, 0, 0);
        end
        
        function rhoxVx(varargin)
            % Three panel figure:
            % Volumetric charge density (rho), Field and potential as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            rho = dfana.calcrho(sol);
            
            figure(211)
            subplot(2, 1, 1)
            dfplot.x2d(sol, x, {rho},{'\rho'},{'-'}, 'Charge density [cm-3]', tarr, xrange, 0, 0);
            
            subplot(2, 1, 2)
            dfplot.x2d(sol, x, {-V},{'V'},{'-'},'-Electrostatic potential [V]', tarr, xrange, 0, 0);
        end
        
        function ELx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            
            figure(22);
            dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'CB', 'VB'},...
                {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)
        end
        
        function ELnpx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);
            
            figure(23)
            subplot(2,1,1);
            dfplot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'CB', 'VB'},...
                {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0);
            
            subplot(2,1,2);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-'}, 'El carrier density [cm-3]', tarr, xrange, 0, 1);
        end
        
        function Vacx(varargin)
            % Potential and ionic charges as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            figure(24)
            subplot(2,1,1);
            dfplot.x2d(sol, x, {V}, {'V'},...
                {'-'}, 'Electro. potential [V]', tarr, xrange, 0, 0);
            
            subplot(2,1,2);
            dfplot.x2d(sol, x, {a, c}, {'a', 'c'},...
                {'-', '-'}, 'Ionic carrier density [cm-3]', tarr, xrange , 0, 0);
        end
        
        function Vnpx(varargin)
            % Potential and ionic charges as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            figure(24)
            subplot(2,1,1);
            dfplot.x2d(sol, x, {V}, {'V'},...
                {'-'}, 'Electro. potential [V]', tarr, xrange, 0, 0);
            
            subplot(2,1,2);
            dfplot.x2d(sol, x, {n, p}, {'n', 'p'},...
                {'-', '-'}, 'Ionic carrier density [cm-3]', tarr, xrange , 0, 1);
        end
        
        function Vionacx(varargin)
            % Electrostatic potential as a function of position
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            Vion = dfana.calcVion(sol);
            Vel = V - Vion;
            
            figure(25)
            subplot(2,1,1);
            dfplot.x2d(sol, x, {V, Vion, Vel}, {'V', 'Vion', 'Vel'},...
                {'--', '.', '-'}, 'Electro. potential [V]', tarr, xrange, 0, 0);
            
            subplot(2,1,2);
            dfplot.x2d(sol, x, {a, c}, {'a', 'c'},...
                {'-', '-'}, 'Ionic carrier density [cm-3]', tarr, xrange , 0, 0);
        end
        
        function Fiont(sol, xpos)
            % Field contribution from ionic charge FION as a function of time at position XPOS
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            
            Fion = dfana.calcFion(sol);
            t = sol.t;
            
            figure(26)
            plot(t, Fion(:,ppos))
            xlabel('Time')
            ylabel('Ion field [Vcm-1]')
        end
        
        function PLx(varargin)
            % Volumetric charge density (rho) as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = dfplot.sortarg(varargin);
            [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);
            
            U = dfana.calcU(sol);
            
            figure(27)
            dfplot.x2d(sol, x, {U.btb}, {'Ubtb'},...
                {'-'}, 'Radiatve recombination rate [cm-3s-1]', tarr, xrange , 0, 0);
        end
        
        function colourblocks(sol, yrange)
            par = sol.par;
            dcum0 = par.dcum0*1e7;   % Convert to nm
            
            % Multicoloured
            % triplets = [1, 0.87, 0.87; 1, 0.9, 0.7; 0.8, 0.9, 1; 1, 0.8, 0.9; 0.8, 1, 0.9; 0.9, 0.8, 1;0.9, 1, 0.8];
            % French flag n-i-p
%             triplets = [0.8, 0.9, 1;...     HTL
%                 1, 1, 1;...                     Active 1
%                 1, 0.87, 0.87; ...                 ETL
%                 0.9, 0.8, 1;...
%                 0.9, 1, 0.8];
            % Mint- 3 layer
            triplets = [0.8, 0.9, 1;...     HTL
                1, 1, 0.7;...                 int1
                1, 1, 1;...                     Active 1
                1, 1, 0.7;...                 int 2
                1, 0.87, 0.87; ...                 ETL
                0.9, 0.8, 1;...
                1, 0.87, 0.87];
            % Homojunction
%             triplets = [0.8, 0.9, 1;...     HTL
%                 1, 1, 0.7;...                 int1
%                 0.9, 1, 0.8;...                     Active 1
%                 1, 1, 0.7;...                 int 2
%                 1, 0.9, 0.7;...                 Active 2
%                 1, 1, 0.7;...                 int 3
%                 1, 0.87, 0.87; ...                 ETL
%                 0.9, 0.8, 1;...
%                 0.9, 1, 0.8];

            % Homojunction - white absorber
%             triplets = [0.8, 0.9, 1;...     HTL
%                 1, 1, 0.7;...                 int1
%                 1, 1, 1;...                     Active 1
%                 1, 1, 1;...                 int 2
%                 1, 1, 1;...                 Active 2
%                 1, 1, 0.7;...                 int 3
%                 1, 0.87, 0.87; ...                 ETL
%                 0.9, 0.8, 1;...
%                 0.9, 1, 0.8];

            % pn
            % triplets = [1, 0.87, 0.87; 0.85, 0.92, 1];
            
            for i =1:length(dcum0)-1
                v = [dcum0(i) yrange(2); dcum0(i+1) yrange(2); dcum0(i+1) yrange(1); dcum0(i) yrange(1)];   % vertices position
                f = [1 2 3 4];    % Faces
                j = i - ((length(triplets)-1)*floor(i/length(triplets)));
                colour = triplets(j,:);
                patch('Faces',f,'Vertices',v,'FaceColor',colour, 'EdgeColor','none');%,'HandleVisibility','off')
            end
            
            hold on
        end
        
        function [sol, tarr, pointtype, xrange] = sortarg(args)
            
            if length(args) == 1
                sol = args{1};
                tarr = sol.t(end);
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(args) == 2
                sol = args{1};
                tarr = args{2};
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(args) == 3
                sol = args{1};
                tarr = args{2};
                xrange = args{3};
                pointtype = 't';
            end
        end
        
        function x2d(sol, xmesh, variables, legstr, linestyle, ylab, tarr, xrange, logx, logy)
            % SOL = solution structure
            % VARIABLES is an array containing the variables for plotting
            % LEGSTR is the legend string
            % YLAB = y-axis label
            % TARR- array of times
            % XRANGE - limits of the plot as a two element vector
            % LOGX, LOGY - switches for log axes
            ax = gca;
            if ishold(ax) == 0
                cla(ax);    % Clear current axis if held
            end
            par = sol.par;
            xnm = xmesh*1e7;
            
            vmin = min(min(cell2mat(variables)));
            vmax = max(max(cell2mat(variables)));
            
            if vmin == 0 && vmax == 0
                vmin = -1;
                vmax = 1;
            end
            vrange = vmax-vmin;
            if isempty(findobj(ax,'Type','patch'))
                switch logy
                    case 0
                        dfplot.colourblocks(sol, [vmin-(vrange), vmax+(vrange)]);
                    case 1
                        dfplot.colourblocks(sol, [0.1*vmin, 10*vmax]);
                end
            end
            
            vmin_tarr = zeros(length(tarr),length(variables));
            vmax_tarr = zeros(length(tarr),length(variables));
            h = zeros(1, length(variables));
            
            for i = 1:length(tarr)
                % find the time
                p1 = find(sol.t <= tarr(i));
                p1 = p1(end);
                for jj = 1:length(variables)
                    vtemp = variables{jj};
                    
                    vmin_tarr(i,jj) = min(vtemp(p1, :));
                    vmax_tarr(i,jj) = max(vtemp(p1, :));
                    
                    h(i,jj) = plot(xnm, variables{jj}(p1, :), char(linestyle(jj)));
                    hold on
                end
            end
            xlabel('Position [nm]')
            ylabel(ylab)
            if logy == 1
                set(gca, 'YScale','log');
            end
            if logx == 1
                set(gca, 'XScale','log');
            end
            if length(variables) == 1
                mystr = [];
                for i = 1:length(tarr)
                    mystr = [mystr, string(['t = ', num2str(tarr(i)), ' s'])];
                end
                lgd = legend(h, mystr);
            else
                lgd = legend(h(1,:), legstr);
            end
            lgd.FontSize = 12;
            xlim([xrange(1), xrange(2)])
            ymin = min(min(vmin_tarr));
            ymax = max(max(vmax_tarr));
            yrange = ymax - ymin;
            if ymin == 0 && ymax == 0
            else
                switch logy
                    case 0
                        if yrange == 0
                            ylim([ymin*0.9, ymax*1.1]);
                        else
                            ylim([ymin-(yrange*0.2), ymax+(yrange*0.2)]);
                        end
                    case 1
                        ylim([0.1*ymin, 10*ymax])
                end
            end
            set(gca, 'Layer', 'top')
            box on
            hold off
        end
    end
end
