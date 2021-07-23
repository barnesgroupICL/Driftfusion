function par = import_properties(par, filepath)
% A function to IMPORT_PROPERTIES from a text file LOCATED at FILEPATH. Each of the listed properties
% is checked to see if it is available in the .CSV file. If it is available, the existing properties
% are overwritten otherwise a warning is displayed. Some entries have
% nested try-ctach statements for backwards compatibility with older
% variable names stored in .csv files. The object properties will still be
% imported with the latest nomenclature.
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Start code
T = readtable(filepath{1,1});   % Reads-in in the external .CSV file to a table T

% Layer type array
try
    par.layer_type = T{:,'layer_type'}';
catch
    warning('No layer type (layer_type) defined in .csv . Using default in PC')
end
% Layer name array
try
    par.stack = T{:,'stack'}';
catch
    warning('No stack (stack) defined in .csv . Using default in PC')
end
% Layer thickness array
try
    par.d = T{:, 'dcell'}';
catch
    try
        par.d = T{:, 'd'}';
    catch
        try
            par.d = T{:, 'thickness'}';
        catch
            warning('No thickness array (thickness) defined in .csv . Using default in PC')
        end
    end
end
% Layer points array
try
    par.layer_points = T{:, 'layer_points'}';
catch
    warning('No layer points array (points) defined in .csv . Using default in PC')
end
% Spatial mesh coefficient for non-linear meshes
try
    par.xmesh_coeff = T{:, 'xmesh_coeff'}';
catch
    warning('No xmesh coefficient array (xmesh_coeff) defined in .csv . Using default in PC')
end
% Electron affinity array
try
    par.EA = T{:, 'EA'}';
catch
    warning('No electron affinity array (EA) defined in .csv . Using default in PC')
end
% Ionisation potential array
try
    par.IP = T{:, 'IP'}';
catch
    warning('No ionisation potential array (IP) defined in .csv . Using default in PC')
end
% Equilibrium Fermi energy array
try
    par.E0 = T{:, 'E0'}';
catch
    warning('No equilibrium Fermi level array (E0) defined in .csv . Using default in PC')
end
% Conduction band effective density of states
try
    par.Nc = T{:, 'Nc'}';
catch
    warning('No conduction band eDOS array (Nc) defined in .csv . Using default in PC')
end
% Valence band effective density of states
try
    par.Nv = T{:, 'Nv'}';
catch
    warning('No valence band eDOS array (Nv) defined in .csv . Using default in PC')
end
% Intrinsic anion density
try
    par.Nani = T{:, 'Nani'}';
catch
    warning('No equilibrium anion density array (Nani) defined in .csv . Using default in PC')
end
% Intrinsic cation density
try
    par.Ncat = T{:, 'Ncat'}';
catch
    warning('No equilibrium cation density array (Ncat) defined in .csv . Using default in PC')
end
% Limiting density of anion states
try
    par.amax = T{:, 'amax'}';
catch
    try
        par.amax = T{:, 'DOSani'}';
    catch
        try
            par.amax = T{:, 'Nion'}';
        catch
            warning('No maximum anion density array (DOSani) defined in .csv . Using default in PC')
        end
    end
end
% Limiting density of cation states
try
    par.cmax = T{:, 'cmax'}';
catch
    try
        par.cmax = T{:, 'DOScat'}';
    catch
        warning('No maximum cation density array (DOScat) defined in .csv . Using default in PC')
    end
end
% Electron mobility
try
    par.mue = T{:, 'mue'}';
catch
    warning('No electron mobility (mue) defined in .csv . Using default in PC')
end
% Hole mobility
try
    par.muh = T{:, 'muh'}';
catch
    warning('No hole mobility (muh) defined in .csv . Using default in PC')
end
% Anion mobility
try
    par.muani = T{:, 'muani'}';
catch
    try
        par.muani = T{:, 'muion'}';
    catch
        warning('No anion mobility (muani) defined in .csv . Using default in PC')
    end
end
% Cation mobility
try
    par.mucat = T{:, 'mucat'}';
catch
    warning('No cation mobility (mucat) defined in .csv . Using default in PC')
end
% Relative dielectric constant
try
    par.epp = T{:, 'epp'}';
catch
    warning('No relative dielectric constant (epp) defined in .csv . Using default in PC')
end
% Uniform volumetric generation rate
try
    par.g0 = T{:, 'g0'}';
catch
    try
        par.g0 = T{:, 'G0'}';
    catch
        warning('No uniform generation rate (g0) defined in .csv . Using default in PC')
    end
end
% Band-to-band recombination coefficient
try
    par.B = T{:, 'krad'}';
catch
    try
        par.B = T{:, 'B'}';
    catch
        warning('No radiative recombinaiton coefficient array (B) defined in .csv . Using default in PC')
    end
end
% Electron SRH time constant
try
    par.taun = T{:, 'taun'}';
catch
    try
        par.taun = T{:, 'taun_SRH'}';
    catch
        warning('No SRH electron lifetime array (taun) defined in .csv . Using default in PC')
    end
end
% Hole SRH time constant
try
    par.taup = T{:, 'taup'}';
catch
    try
        par.taup = T{:, 'taup_SRH'}';
    catch
        warning('No SRH hole lifetime array (taup) defined in .csv . Using default in PC')
    end
end

try
    par.sn = T{:,'sn'}';
catch
    if any(strcmp(par.layer_type, 'interface')) || any(strcmp(par.layer_type, 'junction'))
        warning('No sn value defined in .csv . Using default in PC')
    end
end

try
    par.sp = T{:,'sp'}';
catch
    if any(strcmp(par.layer_type, 'interface')) || any(strcmp(par.layer_type, 'junction'))
        warning('No sp value defined in .csv . Using default in PC')
    end
end

% Electron surface recombination velocity/extraction coefficient LHS
try
    par.sn_l = T{1, 'sn_l'}';
catch
    warning('No sn_l defined in .csv . Using default in PC')
end
% Hole surface recombination velocity/extraction coefficient LHS
try
    par.sp_l = T{1, 'sp_l'}';
catch
    warning('No sp_l defined in .csv . Using default in PC')
end
% Electron surface recombination velocity/extraction coefficient RHS
try
    par.sn_r = T{1, 'sn_r'}';
catch
    warning('No sn_r defined in .csv . Using default in PC')
end
% Hole surface recombination velocity/extraction coefficient RHS
try
    par.sp_r = T{1, 'sp_r'}';
catch
    warning('No sp_r defined in .csv . Using default in PC')
end
% Electrode workfunction LHS
try
    par.Phi_left = T{1, 'Phi_left'};
catch
    try
        par.Phi_left = T{1, 'PhiA'};
    catch
        warning('No Phi_left defined in .csv . Using default in PC')
    end
end
% Electrode workfunction RHS
try
    par.Phi_right = T{1, 'Phi_right'};
catch
    try
        par.Phi_right = T{1, 'PhiC'};
    catch
        warning('No Phi_right defined in .csv . Using default in PC')
    end
end
% SRH Trap energy
try
    par.Et = T{:,'Et'}';
catch
    try
        par.Et = T{:,'Et_bulk'}';
    catch
        warning('No trap energy array (Et) defined in .csv . Using default in PC')
    end
end
% Optical model
try
    par.OM = T{1, 'OM'};
catch
    warning('No optical model (OM) specified in .csv Using default in PC')
end
% Spatial mesh type
try
    par.xmesh_type = T{1, 'xmesh_type'};
catch
    warning('No spatial mesh type (xmesh_type) defined in .csv . Using default in PC')
end
% Illumination side
try
    par.side = T{1, 'side'};
catch
    warning('Illumination side (side) undefined in .csv . Using default in PC')
end
% Number of ionic species
try
    par.N_ionic_species = T{1, 'N_ionic_species'};
catch
    warning('No of ionic species (N_ionic_species) undefined in .csv. Using default in PC')
end
% Layer colours
try
    Red = T{:, 'Red'};
    Green = T{:, 'Green'};
    Blue = T{:, 'Blue'};
    par.layer_colour = [Red,Green,Blue];
catch
    % warning('Layer colours (layer_colour) undefined in .csv. Using default in PC')
end
% Recombination zone location
% Illumination side
if any(strcmp(par.layer_type, 'interface')) || any(strcmp(par.layer_type, 'junction'))
    try
        rec_zone_loc_user = T{:, 'rec_zone_loc'}';
    catch
        warning('Recomination zone location (rec_zone_loc) not defined in .csv . Using default in PC')
        rec_zone_loc_user = strings(1, length(par.stack));
        rec_zone_loc_user(1, 1:length(par.stack)) = deal('auto');
    end
        for i = 1:length(par.stack)
            if any(strcmp(rec_zone_loc_user(i), {'L', 'R', 'C'})) == 1
                par.rec_zone_loc(i) = rec_zone_loc_user(i);
            end
        end
end

end