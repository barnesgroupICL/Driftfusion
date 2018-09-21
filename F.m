classdef F
    % distribution function Class - 
    % calculates carrier densities for different distribution functions -
    % should be renamed!
    properties (Constant)
        % These cannot be altered
        
        % Physical constants
        kB = 8.617330350e-5;       % Boltzmann constant [eV K^-1]
        
    end
    
    methods (Static)
        function n = boltzn(Nc, Ec, Efn, T)
            % Boltzmann carrier density
            
            kT = F.kB*T;
            n = Nc*exp((Efn-Ec)/kT);
            
        end
        
        function p = boltzp(Nv, Ev, Efp, T)
            
            kT = F.kB*T;
            p = Nv*exp((Ev-Efp)/kT);
            
        end
        
        function n = fdn(Nc, Ec, Efn, T)
            % Fermi dirac integral for obtaining electron densities
            % Nc = conduction band density of states
            % Ec = conduction band energy
            % Ef = Fermi level
            % T = temperature
            % See Schubert 2015, pp. 130
            
            for i=1:length(Nc)
                
                kT = F.kB*T;
                
                f = @(E) ((E/kT).^0.5)./(1 + exp((E-Efn(i)+Ec(i))/kT));
                
                n(i) = ((2*Nc(i))/(kT*pi^0.5))*integral(f, 0, Inf);
                
            end
            
        end
        
        function p = fdp(Nv, Ev, Efp, T)
            % Fermi dirac integral for obtaining electron densities
            % Nc = conduction band density of states
            % Ec = conduction band energy
            % Ef = Fermi level
            % T = temperature
            % See Schubert 2015, pp. 130
            
            kT = F.kB*T;
            
            for i=1:length(Nv)
                
                f = @(E) ((E/kT).^0.5).*(1 - 1./(1 + exp((E-Efp(i)+Ev(i))/kT)));
                
                p(i) = imag(((2*Nv(i))/(kT*pi^0.5))*integral(f, -Inf, 0));
                % Not sure why this is the imaginary component- must be something to do with the integral function
            end
            
        end
        
        
    end
    
end

