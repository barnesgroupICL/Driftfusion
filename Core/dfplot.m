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
        
        function EL
            
            
            % Band Diagram
            FH1 = figure(1);
            %set(FigHandle, 'units','normalized','position',[.1 .1 .4 .4]);
            PH1 = subplot(3,1,1);
            %figure(70)
            plot (xnm, Efn(pparr(i),:), '--', xnm, Efp(pparr(i),:), '--', xnm, Ecb(pparr(i), :), xnm, Evb(pparr(i) ,:));
            legend('E_{fn}', 'E_{fp}', 'CB', 'VB');
            set(legend,'FontSize',12);
            xlabel('Position [nm]');
            ylabel('Energy [eV]');
            xlim([xrange(1), xrange(2)]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            grid off;
            drawnow;
            
            hold on
            
            % Final Charge Densities
            %figure(2)
            PH2 = subplot(3,1,2);
            %figure(71)
            semilogy(xnm, n(pparr(i), :), xnm, p(pparr(i), :));
            ylabel('{\itn, p} [cm^{-3}]')
            legend('\itn', '\itp')
            xlabel('Position [nm]')
            xlim([xrange(1), xrange(2)]);
            ylim([1e0, 1e20]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            grid off
            
            hold on
            
            PH3 = subplot(3,1,3);
            %figure(72)
            plot(xnm, (rhoa(pparr(i),:))/1e18, 'black');
            ylabel('{\it\rho a} [x10^{18} cm^{-3}]');
            xlabel('Position [nm]');
            xlim([xrange(1), xrange(2)]);
            %ylim([0, 1.1*(max(sol(pparr(i),:,3))/1e19)]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            grid off
        
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
