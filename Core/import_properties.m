function par = import_properties(par, filepath)
            % A function to overwrite the properties in par with those imported from a
            % text file located in filepath
            T = readtable(filepath{1,1});

            par.layer_type = T{:,'layer_type'}';
            par.stack = T{:,'stack'}';
            par.dcell = T{:, 'thickness'}';
            par.pcell = T{:, 'points'}';
            par.EA = T{:, 'EA'}';
            par.IP = T{:, 'IP'}';
            par.E0 = T{:, 'E0'}';

            par.Nc = T{:, 'Nc'}';
            par.Nv = T{:, 'Nv'}';
            try
                par.Nani = T{:, 'Nion'}';   % Backward compatibility
            end
            try
                par.Nani = T{:, 'Nani'}';
            end
            par.Ncat = T{:, 'Ncat'}';
            try
                par.DOSani = T{:, 'DOSion'}';   % backwards compatibility
            end
            try
                par.DOSani = T{:, 'DOSani'}';
            end
            par.DOScat = T{:, 'DOScat'}';
            par.mue = T{:, 'mue'}';
            par.muh = T{:, 'muh'}';
            try
                par.muani = T{:, 'muion'}';
            end

            try
                par.muani = T{:, 'muani'}';
            end
            par.mucat = T{:, 'mucat'}';
            par.epp = T{:, 'epp'}';

            try
                par.g0 = T{:, 'G0'}';   % For backwards compatibility
            end

            try
                par.g0 = T{:, 'g0'}';
            end

            par.krad = T{:, 'krad'}';
            par.taun = T{:, 'taun_SRH'}';
            par.taup = T{:, 'taup_SRH'}';
            par.sn_l = T{1, 'sn_l'}';
            par.sp_l = T{1, 'sp_l'}';
            par.sn_r = T{1, 'sn_r'}';
            par.sp_r = T{1, 'sp_r'}';

            try
                par.Phi_left = T{1, 'PhiA'};
            end

            try
                par.Phi_left = T{1, 'Phi_left'};    % For backwards compatibility
            end

            try
                par.Phi_right = T{1, 'PhiC'};       % For backwards compatibility
            end

            try
                par.Phi_right = T{1, 'Phi_right'};
            end
            
            try
                par.Rs = T{1, 'Rs'};
            end
            
            try
                par.Et = T{:,'Et'}';
            end

            try
                par.Et = T{:,'Et_bulk'}';   % For backwards compatibility
            end

            try
                par.OM = T{1, 'OM'};
            catch
                warning('No optical model (OM) specified in .csv Using default in PC')
            end

            try
                par.xmesh_type = T{1, 'xmesh_type'};
            catch
                warning('No spatial mesh type (xmesh_type) defined in .csv . Using default in PC')
            end

            try
                par.side = T{1, 'side'};
            catch
                warning('Illumination side (side) undefined in .csv . Using default in PC')
            end

            try
                par.N_ionic_species = T{1, 'N_ionic_species'};
            catch
                warning('No of ionic species (N_ionic_species) undefined in .csv. Using default in PC')
            end

end