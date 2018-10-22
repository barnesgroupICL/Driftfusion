classdef F
    % distribution function Class -
    % calculates carrier densities for different distribution functions -
    % should be renamed!
    properties (Constant)
        % These cannot be altered
        
        % Physical constants
        kB = 8.617330350e-5;       % Boltzmann constant [eV K^-1]
        uplimit = 10;              % Upper limit for the integration [eV]
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
                
                n(i) = real(((2*Nc(i))/(kT*pi^0.5))*integral(fn, 0, F.uplimit));
                
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
                
                p(i) = real(((2*Nv(i))/(kT*pi^0.5))*integral(fp, 0, F.uplimit));
                
            end
            
        end
        
        function Dnfd = Dn_fd_fun(Nc, Ec, Efn, mue, T)
            % Calculates Fermi Dirac diffusion coefficient function
            % Currently uses upper limit for the integral of CB +16eV- using
            % higher values causes problems resulting in Nan results for Dn
            % - requires further investigation
            % Nc = conduction band density of states
            % Ec = conduction band energy
            % Ef = Fermi level
            % T = temperature
            % See Schubert 2015, pp. 130
            
            kB = 8.617330350e-5;       % Boltzmann constant [eV K^-1]
            kT = kB*T;
            e = 1.61917e-19;         % Elementary charge in Coulombs.
            
            for i = 1:length(Efn)
                
                f = @(E) 1./(1 + exp((E-Efn(i)+Ec)/kT));    % Fermi- Dirac function
                dfdE = @(E) exp((E-Efn(i)+Ec)/kT)./(kT*(exp((E-Efn(i)+Ec)/kT)+1).^2);
                
                g = @(E) (E/kT).^0.5;                        % DOS function based on 3D semiconductor
                h = @(E) g(E).*f(E);
                k = @(E) g(E).*dfdE(E);
                
                n(i) = real(((2*Nc)/(kT*pi^0.5))*integral(f, 0, F.uplimit));
                dndE(i) = ((2*Nc)/(kT*pi^0.5))*integral(dfdE, 0, F.uplimit);
            end
            
            Dnfd.Dnfun = mue*(n./dndE);
            Dnfd.n_fd = n;
            Dnfd.Efn = Efn;
        end
        
        function Dpfd = Dp_fd_fun(Nv, Ev, Efp, muh, T)
            % Calculates Fermi Dirac diffusion coefficient function
            % Currently uses upper limit for the integral of CB +16eV- using
            % higher values causes problems resulting in Nan results for Dn
            % - requires further investigation
            % Nc = conduction band density of states
            % Ec = conduction band energy
            % Ef = Fermi level
            % T = temperature
            % See Schubert 2015, pp. 130
            
            % Reflecting the energy makes the integral easier for some
            % reason- doesn't seem to liek integrating from negative
            % infinitiy...
            Efp_trans = 2*Ev-Efp;
            
            kB = 8.617330350e-5;       % Boltzmann constant [eV K^-1]
            kT = kB*T;
            e = 1.61917e-19;         % Elementary charge in Coulombs.
            
            for i = 1:length(Efp)
                
                f = @(E) 1./(1 + exp((E-Efp_trans(i)+Ev)/kT));    % Fermi- Dirac function
                dfdE = @(E) exp((E-Efp_trans(i)+Ev)/kT)./(kT*(exp((E-Efp_trans(i)+Ev)/kT)+1).^2);
                
                g = @(E) (E/kT).^0.5;                        % DOS function based on 3D semiconductor
                h = @(E) g(E).*f(E);
                k = @(E) g(E).*dfdE(E);
                
                p(i) = real(((2*Nv)/(kT*pi^0.5))*integral(f, 0, F.uplimit));
                dpdE(i) = ((2*Nv)/(kT*pi^0.5))*integral(dfdE, 0, F.uplimit);
            end
            
            Dpfd.Dpfun = muh*(p./dpdE);
            Dpfd.p_fd = p;
            Dpfd.Efp = Efp;
        end
        
        function Dsol = D(n, Dnfun, n_fd)
            
            if n > max(n_fd)
                warning('Carrier density is outside of specified limits for the Fermi Dirac integral: Try increasing the Fermi energy limit in PC.BUILDDEV.')
                
            end
            
            % To avoid problems with negative carrier densities
            if n < 0
                n = 0;
            end
            
            if n>min(n_fd)
                pp = find(n_fd <= n);
                pp = pp(end);
                Dsol = Dnfun(1,pp);
            else
                Dsol = min(Dnfun);
            end
            
            %Dsol = interp1(n_fd, Dnfun, n, 'nearest', 'extrap');
            
        end
        
        function Efn_fd = Efn_fd_fun(n, Efn, n_fd)
            
            if n>min(n_fd)
                pp = find(n_fd <= n);
                pp = pp(end);
                Efn_fd = Efn(1,pp);
            else
                Efn_fd = min(Efn);
            end
        end
        
        function Efp_fd = Efp_fd_fun(p, Efp, p_fd)
            
            if p>min(p_fd)
                pp = find(p_fd <= p);
                pp = pp(1);
                Efp_fd = Efp(1,pp);
            else
                Efp_fd = min(Efp);
            end
        end
        
    end
end
