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
            n = zeros(length(Nc));
            
            for i=1:length(Nc)
                
                kT = F.kB*T;
                
                fn = @(E) ((E/kT).^0.5)./(1 + exp((E-Efn(i)+Ec(i))/kT));
                
                n(i) = real(((2*Nc(i))/(kT*pi^0.5))*integral(fn, Ec(i), Inf));
                
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
            p = zeros(length(Nv));
            
            % Reflecting the energy makes the integral easier for some
            % reason- doesn't seem to liek integrating from negative
            % infinitiy...
            Efp = Ev-(Efp-Ev);
            
            for i=1:length(Nv)
                
                fp = @(E) ((E/kT).^0.5)./(1 + exp((E-Efp(i)+Ev(i))/kT));
                
                p(i) = real(((2*Nv(i))/(kT*pi^0.5))*integral(fp, Ev(i), Inf));
                
            end
            
            
            %             for i=1:length(Nv)
            %
            %                 fp = @(E) ((E/kT).^0.5).*(1 - 1./(1 + exp((E-Efp(i)+Ev(i))/kT)));
            %
            %                 p(i) = imag(((2*Nv(i))/(kT*pi^0.5))*integral(fp, -Inf, 0));%Ev(i)));
            %                 % Not sure why this is the imaginary component- must be something to do with the integral function
            %             end
            
        end
        
        function [n,Dn] = Ddfn_numeric(Nc, Ec, Efn, mu, T)
            % Test function for Fermi Dirac diffusion coefficient
            % Curretn uses upper limit for the integral of CB +10 eV- using
            % higher values causes problems resulting in Nan results for Dn
            % Nc = conduction band density of states
            % Ec = conduction band energy
            % Ef = Fermi level
            % T = temperature
            % See Schubert 2015, pp. 130
            
            kB = 8.617330350e-5;       % Boltzmann constant [eV K^-1]
            kT = kB*T;
            e = 1.61917e-19;         % Elementary charge in Coulombs.
            
            for i = 1:length(Efn)
                
                E = Ec:0.01:(10+Ec);
                f = 1./(1 + exp((E-Efn(i)+Ec)/kT));    % Fermi- Dirac function
                dfdE = ((exp((E-Efn(i)+Ec)/kT))./(kT*((exp((E-Efn(i)+Ec)/kT)+1).^2)));
                
                g = (E/kT).^0.5;                        % DOS function based on 3D semiconductor
                h = g.*f;
                k = g.*dfdE;
                
                n(i) = real(((2*Nc)/(kT*pi^0.5))*trapz(E, h));
                dndE(i) = real(((2*Nc)/(kT*pi^0.5))*trapz(E, k));
            end
            
            Dn = mu*(n./dndE);
            
            
        end
        
    end
end
