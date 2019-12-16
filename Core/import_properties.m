function par = import_properties(par, filepath)
% A function to IMPORT_PROPERTIES from a text file LOCATED at FILEPATH. Each of the listed properties
% is checked to see if it is available in the .CSV file. If it is available, the existing properties
% are overwritten otherwise a warning is displayed.

T = readtable(filepath{1,1});   % Reads-in in the external .CSV file to a table T

try
    par.layer_type = T{:,'layer_type'}';
catch
    warning('No layer type (layer_type) defined in .csv . Using default in PC')
end
try
    par.stack = T{:,'stack'}';
catch
    warning('No stack (stack) defined in .csv . Using default in PC')
end
try
    par.dcell = T{:, 'thickness'}';
catch
    warning('No thickness array (thickness) defined in .csv . Using default in PC')
end
try
    par.layer_points = T{:, 'layer_points'}';
catch
    warning('No layer points array (points) defined in .csv . Using default in PC')
end
try
    par.EA = T{:, 'EA'}';
catch
    warning('No electron affinity array (EA) defined in .csv . Using default in PC')
end
try
    par.IP = T{:, 'IP'}';
catch
    warning('No ionisation potential array (IP) defined in .csv . Using default in PC')
end
try
    par.E0 = T{:, 'E0'}';
catch
    warning('No equilibrium Fermi level array (E0) defined in .csv . Using default in PC')
end
try
    par.Nc = T{:, 'Nc'}';
catch
    warning('No conduction band eDOS array (Nc) defined in .csv . Using default in PC')
end
try
    par.Nv = T{:, 'Nv'}';
catch
    warning('No valence band eDOS array (Nv) defined in .csv . Using default in PC')
end
try
    par.Nani = T{:, 'Nani'}';
catch
    warning('No equilibrium anion density array (Nani) defined in .csv . Using default in PC')
end
try
    par.Ncat = T{:, 'Ncat'}';
catch
    warning('No equilibrium cation density array (Ncat) defined in .csv . Using default in PC')
end
try
    par.DOSani = T{:, 'DOSani'}';
catch
    warning('No maximum anion density array (DOSani) defined in .csv . Using default in PC')
end
try
    par.DOScat = T{:, 'DOScat'}';
catch
    warning('No maximum cation density array (DOScat) defined in .csv . Using default in PC')
end
try
    par.mue = T{:, 'mue'}';
catch
    warning('No electron mobility (mue) defined in .csv . Using default in PC')
end
try
    par.muh = T{:, 'muh'}';
catch
    warning('No hole mobility (muh) defined in .csv . Using default in PC')
end

try
    par.muani = T{:, 'muani'}';
catch
    warning('No anion mobility (muani) defined in .csv . Using default in PC')
end
try
    par.mucat = T{:, 'mucat'}';
catch
    warning('No cation mobility (mucat) defined in .csv . Using default in PC')
end
try
    par.epp = T{:, 'epp'}';
catch
    warning('No relative dielectric constant (epp) defined in .csv . Using default in PC')
end

try
    par.g0 = T{:, 'g0'}';
catch
    warning('No uniform generation rate (g0) defined in .csv . Using default in PC')
end
try
    par.krad = T{:, 'krad'}';
catch
    warning('No radiative recombinaiton coefficient array (krad) defined in .csv . Using default in PC')
end

try
    par.taun = T{:, 'taun'}';
catch
    warning('No SRH electron lifetime array (taun) defined in .csv . Using default in PC')
end
try
    par.taup = T{:, 'taup'}';
catch
    warning('No SRH hole lifetime array (taup) defined in .csv . Using default in PC')
end
try
    par.sn_l = T{1, 'sn_l'}';
catch
    warning('No sn_l defined in .csv . Using default in PC')
end
try
    par.sp_l = T{1, 'sp_l'}';
catch
    warning('No sp_l defined in .csv . Using default in PC')
end
try
    par.sn_r = T{1, 'sn_r'}';
catch
    warning('No sn_r defined in .csv . Using default in PC')
end
try
    par.sp_r = T{1, 'sp_r'}';
catch
    warning('No sp_r defined in .csv . Using default in PC')
end

try
    par.Phi_left = T{1, 'Phi_left'};
catch
    warning('No Phi_left defined in .csv . Using default in PC')
end
try
    par.Phi_right = T{1, 'Phi_right'};
catch
    warning('No Phi_right defined in .csv . Using default in PC')
end
try
    par.Et = T{:,'Et'}';
catch
    warning('No trap energy array (Et) defined in .csv . Using default in PC')
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

%% Additional checks for backwards compatibility- old variable names
% No warning required otherwise will constantly warn when building
% parameters object
try
    par.Nani = T{:, 'Nion'}';
end
try
    par.DOSani = T{:, 'DOSion'}';
end
try
    par.muani = T{:, 'muion'}';
end
try
    par.Et = T{:,'Et_bulk'}';
end
try
    par.taun = T{:, 'taun_SRH'}';
end
try
    par.taup = T{:, 'taup_SRH'}';
end
try
    par.g0 = T{:, 'G0'}';
end
try
    par.Phi_left = T{1, 'PhiA'};
end
try
    par.Phi_right = T{1, 'PhiC'}; 
end
try
    par.layer_points = T{:, 'points'}';
end
end