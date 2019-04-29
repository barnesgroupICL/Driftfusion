classdef dfana_class
    
    methods (Static)
        
        function [Efn, Efp] = QFLs(sol)
            % u is the solution structure
            % Simple structure names
            u = sol.u;
            par = sol.par;
            x = sol.x;
            t = sol.t;
            dev = par.dev;
                        
            %% ANALYSIS %%
            xnm = x*1e7;    % x in nm for plotting
            
            % split the solution into its component parts (e.g. electrons, holes and efield)
            n = u(:,:,1);
            p = u(:,:,2);
            a = u(:,:,3);
            V = u(:,:,4);
                      
            % Set minority carrier densities to zero to avoid errors
            nmod = n;
            pmod = p;
            
            % Create 2D matrices for multiplication with solutions
            EAmat = repmat(par.dev.EA, length(t), 1);
            IPmat = repmat(par.dev.IP, length(t), 1);
            muemat = repmat(par.dev.mue, length(t), 1);
            muhmat = repmat(par.dev.muh, length(t), 1);
            muionmat = repmat(par.dev.muion, length(t), 1);
            NAmat = repmat(par.dev.NA, length(t), 1);
            NDmat = repmat(par.dev.ND, length(t), 1);
            Ncmat = repmat(par.dev.Nc, length(t), 1);
            Nvmat = repmat(par.dev.Nv, length(t), 1);
            Nionmat = repmat(par.dev.Nion, length(t), 1);
            eppmat = repmat(par.dev.epp, length(t), 1);
            nimat = repmat(par.dev.ni, length(t), 1);
            kradmat = repmat(par.dev.krad, length(t), 1);
            taunmat = repmat(par.dev.taun, length(t), 1);
            taupmat = repmat(par.dev.taup, length(t), 1);
            ntmat = repmat(par.dev.nt, length(t), 1);
            ptmat = repmat(par.dev.pt, length(t), 1);
            
            Ecb = EAmat-V;                                 % Conduction band potential
            Evb = IPmat-V;                                 % Valence band potential
            
            Efn = zeros(size(n,1), size(n,2));
            Efp = zeros(size(n,1), size(n,2));
                       
            if par.stats == 'Fermi'
                
                for i = 1:size(n,1)           % time
                    for j = 1:size(n,2)       % position
                        Efn(i,j) = F.Efn_fd_fun(nmod(i,j), par.dev.Efn(j,:),  par.dev.n_fd(j,:));
                        Efp(i,j) = F.Efp_fd_fun(pmod(i,j), par.dev.Efp(j,:),  par.dev.p_fd(j,:));
                    end
                end
                Efn = Efn-V;
                Efp = Efp-V;
                
            elseif par.stats == 'Boltz'
                Efn = real(Ecb+(par.kB*par.T/par.q)*log(nmod./Ncmat));        % Electron quasi-Fermi level
                Efp = real(Evb-(par.kB*par.T/par.q)*log(pmod./Nvmat));        % Hole quasi-Fermi level
            end
            
        end
        
    end
    
end
