classdef dfplot
    
    methods (Static)
        
        function JV(JV, option)
            % JV - a soultion from doJV
            % OPTION - 1 = dark only, 2 = light only, 3 = dark & light
            % JV is a structure containing dark and illuminated JVs
            
            if option == 1 || option == 3
                JV.dk.f.par.JV = 1;
                JV.dk.r.par.JV = 1;
                
                dfana(JV.dk.f);
                figure(11)
                %xlim([-0.2, 1.15])
                ylim([-30e-3, 10e-3]);
                hold on
                dfana(JV.dk.r);
                
                JV.dk.f.par.JV = 0;
                JV.dk.r.par.JV = 0;
            end
            
            if option == 2 || option == 3
                JV.ill.f.par.JV = 1;
                JV.ill.r.par.JV = 1;
                
                dfana(JV.ill.f);
                
                figure(11)
                %xlim([-0.2, 1.15])
                ylim([-30e-3, 10e-3]);
                hold on
                dfana(JV.ill.r);
                
                JV.ill.f.par.JV = 0;
                JV.ill.r.par.JV = 0;
            end
            
            hold off
            
        end
        
        function EL(varargin)
            % Energy Level diagram plotter
            % SOL = the solution structure
            % TIME = the time that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            
            % tarr is a time time array for the time you wish to plot
            if length(varargin) == 1
                sol = varargin{1};
                time = sol.t(end);
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 2
                sol = varargin{1};
                time = varargin{2};
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 3
                sol = varargin{1};
                time = varargin{2};
                xrange = varargin{3};
                pointtype = 't';
            end
            
            % find the time
            p1 = find(sol.t <= time);
            p1 = p1(end);
            
            % Call dfana_class to obtain band energies and QFLs
            [u,t,x,par,dev,n,p,a,V] = dfana_class.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana_class.QFLs(sol);
            
            % calc static background ion matrix
            Nionmat = repmat(par.dev.Nion, length(t), 1);
            
            xnm = x*1e7;    % x in nm for plotting
            
            % Band Diagram
            FH1 = figure(1);
            PH1 = subplot(3,1,1);
            plot (xnm, Efn(p1,:), '--', xnm, Efp(p1,:), '--', xnm, Ecb(p1, :), xnm, Evb(p1 ,:));
            legend('E_{fn}', 'E_{fp}', 'CB', 'VB');
            set(legend,'FontSize',12);
            xlabel('Position [nm]');
            ylabel('Energy [eV]');
            xlim([xrange(1), xrange(2)]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            grid off;
            drawnow;
            
            % Final Charge Densities
            PH2 = subplot(3,1,2);
            semilogy(xnm, n(p1, :), xnm, p(p1, :));
            ylabel('{\itn, p} [cm^{-3}]')
            legend('\itn', '\itp')
            xlabel('Position [nm]')
            xlim([xrange(1), xrange(2)]);
            ylim([1e0, 1e20]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            grid off
            
            % Ionic space charge density
            PH3 = subplot(3,1,3);
            plot(xnm, a(p1,:)-Nionmat, 'black');
            ylabel('{\it\rho a} [cm^{-3}]');
            xlabel('Position [nm]');
            xlim([xrange(1), xrange(2)]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            grid off
            
        end
        
        function Jt(sol)
            % Currents as a function of time
            t = sol.t;
            [j, J] = dfana_class.calcJ(sol);
            
            figure(2);
            plot(t, J.n(:, end),t, J.p(:, end),t, J.a(:, end),t, J.tot(:, end));
            legend('Jn', 'Jp', 'Ja', 'Jtotal')
            xlabel('time [s]');
            ylabel('J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
            grid off;
            drawnow;
        end
        
        function Jx(varargin)
            % Plots the currents
            % SOL = the solution structure
            % TIME = the time that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            % tarr is a time time array for the time you wish to plot
            if length(varargin) == 1
                sol = varargin{1};
                time = sol.t(end);
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 2
                sol = varargin{1};
                time = varargin{2};
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 3
                sol = varargin{1};
                time = varargin{2};
                xrange = varargin{3};
                pointtype = 't';
            end
            
            % find the time
            p1 = find(sol.t <= time);
            p1 = p1(end);
            
            xnm = sol.x*1e7;
            [j, J] = dfana_class.calcJ(sol);
            
            % electron and hole currents as function of position from continuity
            figure(3)
            plot(xnm, J.n(p1, :), xnm, J.p(p1, :), xnm, J.a(p1, :), xnm, J.disp(p1, :), xnm, J.tot(p1, :))
            xlabel('Position [nm]')
            ylabel('J [A]')
            legend('Jn', 'Jp', 'Ja', 'Jdisp', 'Jtot')
            xlim([xrange(1), xrange(2)])
            
        end
        
        
        % multiplot 1
        function mp1(varargin)
            
            % tarr is a time time array for the time you wish to plot
            if length(varargin) == 1
                solstruct = varargin{1};
                tarr = solstruct.t(end);
                pointtype = 't';
            elseif length(varargin) == 2
                solstruct = varargin{1};
                tarr = varargin{2};
                pointtype = 't';
            elseif length(varargin) == 3
                solstruct = varargin{1};
                pointtype = varargin{2};
                tarr = varargin{3};
            end
            
        end
        
    end
    
end
