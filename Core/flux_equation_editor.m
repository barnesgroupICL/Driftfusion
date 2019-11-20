function flux = flux_equation_editor(u, dudx, Dn, Dp, i, par)

%% Constants
kB = par.kB;
T = par.T;

%% Variables
n = u(1);
p = u(2);
c = u(3);
V = u(4);
if par.N_ionic_species == 2
    a = u(5);
end

%% Gradients
dndx = dudx(1);
dpdx = dudx(2);
dcdx = dudx(3);
dVdx = dudx(4);
if par.N_ionic_species == 2
    dadx = dudx(5);
end

%% Switches and accelerator coefficients
mobset = par.mobset;        % Electronic carrier transport switch
mobseti = par.mobseti;      % Ionic carrier transport switch
K_cation = par.K_cation;    % Cation transport rate multiplier
K_anion = par.K_anion;      % Anion transport rate multiplier

%% Device parameters
device = par.dev_ihalf;
mue = device.mue(i);      % Electron mobility
muh = device.muh(i);      % Hole mobility
mucat = device.mucat(i);  % Cation mobility
muani = device.muani(i);  % Anion mobility
Nc = device.Nc(i);        % Conduction band effective density of states
Nv = device.Nv(i);        % Valence band effective density of states
Ncat = device.DOScat(i);  % Cation density upper limit
Nani = device.DOSani(i);  % Anion density upper limit
gradNc = device.Nc(i);    % Conduction band effective density of states gradient
gradNv = device.Nv(i);    % Valence band effective density of states gradient
gradEA = device.EA(i);    % Electron Affinity gradient
gradIP = device.IP(i);    % Ionisation Potential gradient
epp = device.epp(i);      % Dielectric constant

%% Electron flux
jn = mobset*(mue*n*(-dVdx+gradEA) + (Dn*(dndx-((n/Nc)*gradNc))));

%% Hole flux
jp = mobset*(muh*p*(dVdx-gradIP) + (Dp*(dpdx-((p/Nv)*gradNv))));

%% Cation flux
jc = K_cation*mobseti*(mucat*(c*dVdx + kB*T*(dcdx+(c*(dcdx/(Ncat-c))))));       % Nerst-Planck-Poisson approach ref: Borukhov 1997

%% Electrostatic potential 'flux'
jV = (epp/max(par.epp))*dVdx;

if par.N_ionic_species == 2
    ja = K_anion*mobseti*(muani*(a*-dVdx + kB*T*(dadx+(a*(dadx/(Nani-a))))));
end

switch par.N_ionic_species
    case 1
        flux = [jn; jp; jc; jV];
    case 2
         flux = [jn; jp; jc; jV; ja];
end

end

