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
    layer_type = T{:,'layer_type'}';
    if strcmp(layer_type{1}, 'electrode')
        start_row = 2;
        end_row = length(layer_type) - 1;
        par.layer_type = layer_type(start_row:end_row);
    else
        start_row = 1;
        end_row = length(layer_type);
    end
    par.layer_type = layer_type(start_row:end_row);
catch
    error('No layer type (layer_type) defined in .csv . layer_type must be defined when using .csv input file')
end
% Layer name array
try
    material = T{:,'material'}';
    par.material = material(start_row:end_row);
catch
    try
        material = T{:,'stack'}';
        par.material = material(start_row:end_row);
    catch
        warning('No material array (material) defined in .csv . Using default in PC')
    end
end
% Layer thickness array
try
    d = T{:, 'dcell'}';
    par.d = d(start_row:end_row);
catch
    try
        d = T{:, 'd'}';
        par.d = d(start_row:end_row);
    catch
        try
            d = T{:, 'thickness'}';
            par.d = d(start_row:end_row);
        catch
            warning('No thickness array (thickness) defined in .csv . Using default in PC')
        end
    end
end
% Layer points array
try
    layer_points = T{:, 'layer_points'}';
    par.layer_points = layer_points(start_row:end_row);
catch
    warning('No layer points array (points) defined in .csv . Using default in PC')
end
% Spatial mesh coefficient for non-linear meshes
try
    xmesh_coeff = T{:, 'xmesh_coeff'}';
    par.xmesh_coeff = xmesh_coeff(start_row:end_row);
catch
    warning('No xmesh coefficient array (xmesh_coeff) defined in .csv . Using default in PC')
end
% Electron affinity array
try
    Phi_EA = T{:, 'Phi_EA'}';
    par.Phi_EA = Phi_EA(start_row:end_row);
catch
    try
        Phi_EA = T{:, 'EA'}';
        par.Phi_EA = Phi_EA(start_row:end_row);
    catch
        warning('No electron affinity array (Phi_EA) defined in .csv . Using default in PC')
    end
end
% Ionisation potential array
try
    Phi_IP = T{:, 'Phi_IP'}';
    par.Phi_IP = Phi_IP(start_row:end_row);
catch
    try
        Phi_IP = T{:, 'IP'}';
        par.Phi_IP = Phi_IP(start_row:end_row);
    catch
        warning('No ionisation potential array (Phi_IP) defined in .csv . Using default in PC')
    end
end
% Equilibrium Fermi energy array
try
    EF0 = T{:, 'EF0'}';
    par.EF0 = EF0(start_row:end_row);
    if strcmp(layer_type{1}, 'electrode')
        par.Phi_left = EF0(1);
        par.Phi_right = EF0(end);
    end
catch
    try
        EF0 = T{:, 'E0'}';
        par.EF0 = EF0(start_row:end_row);
        if strcmp(layer_type{1}, 'electrode')
            par.Phi_left = EF0(1);
            par.Phi_right = EF0(end);
        end
    catch
        warning('No equilibrium Fermi level array (EF0) defined in .csv . Using default in PC')
    end
end
% Conduction band effective density of states
try
    Nc = T{:, 'Nc'}';
    par.Nc = Nc(start_row:end_row);
catch
    warning('No conduction band eDOS array (Nc) defined in .csv . Using default in PC')
end
% Valence band effective density of states
try
    Nv = T{:, 'Nv'}';
    par.Nv = Nv(start_row:end_row);
catch
    warning('No valence band eDOS array (Nv) defined in .csv . Using default in PC')
end
% Intrinsic anion density
try
    Nani = T{:, 'Nani'}';
    par.Nani = Nani(start_row:end_row);
catch
    warning('No equilibrium anion density array (Nani) defined in .csv . Using default in PC')
end
% Intrinsic cation density
try
    Ncat = T{:, 'Ncat'}';
    par.Ncat = Ncat(start_row:end_row);
catch
    try
        Ncat = T{:, 'Nion'}';
        par.Ncat = Ncat(start_row:end_row);
    catch
        warning('No equilibrium cation density array (Ncat) defined in .csv . Using default in PC')
    end
end
% Limiting density of anion states
try
    a_max = T{:, 'a_max'}';
    par.a_max = a_max(start_row:end_row);
catch
    try
        a_max = T{:, 'DOSani'}';
        par.a_max = a_max(start_row:end_row);
    catch
        try
            a_max = T{:, 'amax'}';
            par.a_max = a_max(start_row:end_row);
        catch
            warning('No maximum anion density array (a_max) defined in .csv . Using default in PC')
        end
    end
end
% Limiting density of cation states
try
    c_max = T{:, 'c_max'}';
    par.c_max = c_max(start_row:end_row);
catch
    try
        c_max = T{:, 'DOScat'}';
        par.c_max = c_max(start_row:end_row);
    catch
        warning('No maximum cation density array (c_max) defined in .csv . Using default in PC')
    end
end
% Electron mobility
try
    mu_n = T{:, 'mu_n'}';
    par.mu_n = mu_n(start_row:end_row);
catch
    try
        mu_n = T{:, 'mue'}';
        par.mu_n = mu_n(start_row:end_row);
    catch
        warning('No electron mobility (mu_n) defined in .csv . Using default in PC')
    end
end
% Hole mobility
try
    mu_p = T{:, 'mu_p'}';
    par.mu_p = mu_p(start_row:end_row);
catch
    try
        mu_p = T{:, 'muh'}';
        par.mu_p = mu_p(start_row:end_row);
    catch
        warning('No hole mobility (mu_p) defined in .csv . Using default in PC')
    end
end
% Anion mobility
try
    mu_a = T{:, 'mu_a'}';
    par.mu_a = mu_a(start_row:end_row);
catch
    try
        mu_a = T{:, 'muani'}';
        par.mu_a = mu_a(start_row:end_row);
    catch
        warning('No anion mobility (mu_a) defined in .csv . Using default in PC')
    end
end
% Cation mobility
try
    mu_c = T{:, 'mu_c'}';
    par.mu_c = mu_c(start_row:end_row);
catch
    try
        mu_c = T{:, 'mucat'}';
        par.mu_c = mu_c(start_row:end_row);
    catch
        warning('No cation mobility (mu_c) defined in .csv . Using default in PC')
    end
end
% Relative dielectric constant
try
    epp = T{:, 'epp'}';
    par.epp = epp(start_row:end_row);
catch
    warning('No relative dielectric constant (epp) defined in .csv . Using default in PC')
end
% Uniform volumetric generation rate
try
    g0 = T{:, 'g0'}';
    par.g0 = g0(start_row:end_row);
catch
    try
        g0 = T{:, 'G0'}';
        par.g0 = g0(start_row:end_row);
    catch
        warning('No uniform generation rate (g0) defined in .csv . Using default in PC')
    end
end
% Band-to-band recombination coefficient
try
    B = T{:, 'krad'}';
    par.B = B(start_row:end_row);
catch
    try
        B = T{:, 'B'}';
        par.B = B(start_row:end_row);
    catch
        warning('No radiative recombinaiton coefficient array (B) defined in .csv . Using default in PC')
    end
end
% Electron SRH time constant
try
    taun = T{:, 'taun'}';
    par.taun = taun(start_row:end_row);
catch
    try
        taun = T{:, 'taun_SRH'}';
        par.taun = taun(start_row:end_row);
    catch
        warning('No SRH electron lifetime array (taun) defined in .csv . Using default in PC')
    end
end
% Hole SRH time constant
try
    taup = T{:, 'taup'}';
    par.taup = taup(start_row:end_row);
catch
    try
        taup = T{:, 'taup_SRH'}';
        par.taup = taup(start_row:end_row);
    catch
        warning('No SRH hole lifetime array (taup) defined in .csv . Using default in PC')
    end
end

try
    sn = T{:,'sn'}';
    par.sn = sn(start_row:end_row);
    if strcmp(layer_type{1}, 'electrode')
        par.sn_l = sn(1);
        par.sn_r = sn(end);
    end
catch
    if any(strcmp(par.layer_type, 'interface')) || any(strcmp(par.layer_type, 'junction'))
        warning('No sn value defined in .csv . Using default in PC')
    end
end

try
    sp = T{:,'sp'}';
    par.sp = sp(start_row:end_row);
    if strcmp(layer_type{1}, 'electrode')
        par.sp_l = sp(1);
        par.sp_r = sp(end);
    end
catch
    if any(strcmp(par.layer_type, 'interface')) || any(strcmp(par.layer_type, 'junction'))
        warning('No sp value defined in .csv . Using default in PC')
    end
end

if strcmp(layer_type{1}, 'electrode') == 0
    % Electron surface recombination velocity/extraction coefficient LHS
    try
        sn_l = T{1, 'sn_l'}';
        par.sn_l = sn_l;
    catch
        warning('No sn_l defined in .csv . Using default in PC')
    end
    % Hole surface recombination velocity/extraction coefficient LHS
    try
        sp_l = T{1, 'sp_l'}';
        par.sp_l = sp_l;
    catch
        warning('No sp_l defined in .csv . Using default in PC')
    end
    % Electron surface recombination velocity/extraction coefficient RHS
    try
        sn_r = T{1, 'sn_r'}';
        par.sn_r = sn_r;
    catch
        warning('No sn_r defined in .csv . Using default in PC')
    end
    % Hole surface recombination velocity/extraction coefficient RHS
    try
        sp_r = T{1, 'sp_r'}';
        par.sp_r = sp_r;
    catch
        warning('No sp_r defined in .csv . Using default in PC')
    end
    % Electrode workfunction LHS
    try
        Phi_left = T{1, 'Phi_left'};
        par.Phi_left = Phi_left;
    catch
        try
            Phi_left = T{1, 'PhiA'};
            par.Phi_left = Phi_left;
        catch
            warning('No Phi_left defined in .csv . Using default in PC')
        end
    end
    % Electrode workfunction RHS
    try
        Phi_right = T{1, 'Phi_right'};
        par.Phi_right = Phi_right;
    catch
        try
            Phi_right = T{1, 'PhiC'};
            par.Phi_right = Phi_right;
        catch
            warning('No Phi_right defined in .csv . Using default in PC')
        end
    end
end
% SRH Trap energy
try
    Et = T{:,'Et'}';
    par.Et = Et(start_row:end_row);
catch
    try
        Et = T{:,'Et_bulk'}';
        par.Et = Et(start_row:end_row);
    catch
        warning('No trap energy array (Et) defined in .csv . Using default in PC')
    end
end
% Optical model
try
    par.optical_model = T{1, 'optical_model'};
catch
    try
        par.optical_model = T{1, 'OM'};
    catch
        warning('No optical model (optical_model) specified in .csv Using default in PC')
    end
end
% Spatial mesh type
try
    xmesh_type = T{1, 'xmesh_type'};
    if isa(xmesh_type, 'double')
        switch xmesh_type
            case 4
                par.xmesh_type = 'linear';
            case 5
                par.xmesh_type = 'erf-linear';
            otherwise
                error('xmesh_type not recognized')
        end
    else
        par.xmesh_type = xmesh_type{1};
    end
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
    par.layer_colour = [Red(start_row:end_row),Green(start_row:end_row),Blue(start_row:end_row)];
catch
    % warning('Layer colours (layer_colour) undefined in .csv. Using default in PC')
end

% Recombination zone location
if any(strcmp(par.layer_type, 'interface')) || any(strcmp(par.layer_type, 'junction'))
    vsr_zone_loc_user = cell(1, length(par.material));
    par.vsr_zone_loc = cell(1, length(par.material));
    vsr_zone_loc_auto = locate_vsr_zone(par);
    try
        vsr_zone_loc_user = T{:, 'vsr_zone_loc'}';
        vsr_zone_loc_user = vsr_zone_loc_user(start_row:end_row);
    catch
        warning('Recomination zone location (vsr_zone_loc) not defined in .csv . Using auto defined')
        par.vsr_zone_loc = vsr_zone_loc_auto;
    end
        for i = 1:length(par.material)
            if any(strcmp(vsr_zone_loc_user(i), {'L','C','R'})) == 1
                par.vsr_zone_loc(i) = vsr_zone_loc_user(i);
            elseif strcmp(vsr_zone_loc_user(i), {'auto'}) == 1
                par.vsr_zone_loc(i) = vsr_zone_loc_auto(i);
            end
        end
end

end