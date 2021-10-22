function par = import_properties(par, filepath)
% A function to IMPORT_PROPERTIES from a text file LOCATED at FILEPATH. Each of the listed properties
% is checked to see if it is available in the .CSV file. If it is available, the existing properties
% are overwritten otherwise a warning is displayed. Some entries have
% nested try-ctach statements for backwards compatibility with older
% variable names stored in .csv files. The object properties will still be
% imported with the latest nomenclature.
% 
% IMPORT_SINGLE_PROPERTY checks for the existence of column headers in the .CSV
% with POSSIBLE_HEADERS. Multiple names for variables in the .CSV are given
% for backwards compatibility, however once read-in the properties take the
% naming convention of the current version i.e. this is simply to allow old
% parameters files to be read-in.
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
% Material name array
par.material = import_single_property(T, {'material', 'stack'}, start_row, end_row);
% Layer thickness array
par.d = import_single_property(T, {'dcell', 'd', 'thickness'}, start_row, end_row);
% Layer points array
par.layer_points = import_single_property(T, {'layer_points', 'points'}, start_row, end_row);
% Spatial mesh coefficient for non-linear meshes
par.xmesh_coeff = import_single_property(T, {'xmesh_coeff'}, start_row, end_row);
% Electron affinity array
par.Phi_EA = import_single_property(T, {'Phi_EA', 'EA'}, start_row, end_row);
% Ionisation potential array
par.Phi_IP = import_single_property(T, {'Phi_IP', 'IP'}, start_row, end_row);
% Equilibrium Fermi energy array
% SRH Trap energy
par.Et = import_single_property(T, {'Et', 'Et_bulk'}, start_row, end_row);
if strcmp(layer_type{1}, 'electrode')
    EFO = import_single_property(T, {'EF0', 'E0'}, 1, length(layer_type));
    par.EF0 = EFO(start_row:end_row);
    par.Phi_left = EFO(1);
    par.Phi_right = EFO(end);
else
    par.EF0 = import_single_property(T, {'EF0', 'E0'}, start_row, end_row);
end
% Conduction band effective density of states
par.Nc = import_single_property(T, {'Nc', 'Ncb', 'NC', 'NCB'}, start_row, end_row);
% Valence band effective density of states
par.Nv = import_single_property(T, {'Nv', 'Nvb', 'NV', 'NVB'}, start_row, end_row);
% Intrinsic anion density
par.Nani = import_single_property(T, {'Nani'}, start_row, end_row);
% Intrinsic cation density
par.Ncat = import_single_property(T, {'Ncat', 'Nion'}, start_row, end_row);
% Limiting density of anion states
par.a_max = import_single_property(T, {'a_max', 'amax', 'DOSani'}, start_row, end_row);
% Limiting density of cation states
par.c_max = import_single_property(T, {'c_max', 'cmax', 'DOScat'}, start_row, end_row);
% Electron mobility
par.mu_n = import_single_property(T, {'mu_n', 'mun', 'mue', 'mu_e'}, start_row, end_row);
% Hole mobility
par.mu_p = import_single_property(T, {'mu_p', 'mup', 'muh', 'mu_h'}, start_row, end_row);
% Anion mobility
par.mu_a = import_single_property(T, {'mu_a', 'mua',  'mu_ani', 'muani'}, start_row, end_row);
% Cation mobility
par.mu_c = import_single_property(T, {'mu_c', 'muc',  'mu_cat', 'mucat'}, start_row, end_row);
% Relative dielectric constant
par.epp = import_single_property(T, {'epp', 'eppr'}, start_row, end_row);
% Uniform volumetric generation rate
par.g0 = import_single_property(T, {'g0', 'G0'}, start_row, end_row);
% Band-to-band recombination coefficient
par.B = import_single_property(T, {'B', 'krad', 'kbtb'}, start_row, end_row);
% Electron SRH time constant
par.taun = import_single_property(T, {'taun', 'taun_SRH'}, start_row, end_row);
% Hole SRH time constant
par.taup = import_single_property(T, {'taup', 'taup_SRH'}, start_row, end_row);
% Electron and hole surface recombination velocities
if strcmp(layer_type{1}, 'electrode')
    sn = import_single_property(T, {'sn'}, 1, length(layer_type));
    par.sn = sn(start_row:end_row);
    par.sn_l = sn(1);
    par.sn_r = sn(end);
    
    sp = import_single_property(T, {'sp'}, 1, length(layer_type));
    par.sp = sp(start_row:end_row);
    par.sp_l = sp(1);
    par.sp_r = sp(end); 
else
    par.sn = import_single_property(T, {'sn'}, start_row, end_row);
    par.sp = import_single_property(T, {'sp'}, start_row, end_row);
end
   
if strcmp(layer_type{1}, 'electrode') == 0
    % Electron surface recombination velocity/extraction coefficient LHS
    par.sn_l = import_single_property(T, {'sn_l', 'snl'}, 1, 1);
    par.sn_r = import_single_property(T, {'sn_r', 'snr'}, 1, 1);
    % Hole surface recombination velocity/extraction coefficient LHS
    par.sp_l = import_single_property(T, {'sp_l', 'spl'}, 1, 1);
    par.sp_r = import_single_property(T, {'sp_r', 'spr'}, 1, 1);
    % Electrode workfunction LHS
    par.Phi_left = import_single_property(T, {'Phi_left', 'Phi_l', 'PhiA'}, 1, 1);
    % Electrode workfunction RHS
    par.Phi_right = import_single_property(T, {'Phi_right', 'Phi_r', 'PhiC'}, 1, 1);
end

% Optical model
try
    optical_model = T{1, 'optical_model'};
        if isa(optical_model, 'double')
            switch optical_model
                case 0
                   par.optical_model = 'uniform';
                case 1
                   par.optical_model = 'beer-lambert';
            end
        else
            par.optical_model = optical_model{1};
        end
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

    function parameter = import_single_property(T, possible_headers, start_row, end_row)
        error_checker = zeros(1, length(possible_headers));
        
        for i = 1:length(possible_headers)
            try
                temp_param = T{: , possible_headers(i)}';
                parameter = temp_param(start_row:end_row);
            catch
                error_checker(i) = 1;
            end
        end
        
        if all(error_checker)
            warning(['No column headings match', possible_headers, ', using default in PC.'])
        end
        

    end
end