function fromSteadyStateStructToTxt(struct, prefix)
% save the main data from a stabilized solution created by pindrift to
% txt files, ideally easy to import with Origin (from OriginLab)
% example:
%   fromSteadyStateStructToTxt(ssol_i_1S_SR, 'pinParams_ssol_i_1S_SR')

p = struct.p;

%% create header

headerIonic = ["x", "a"];
headerBand = ["x", "CB", "Efn", "Efp", "VB"];

%% get measure units

unitsIonic = ["nm", "cm\+(-3)"];
unitsBand = ["nm", "eV", "eV", "eV", "eV"];

%% get data

x = struct.x;
xnm = x * 1e7;

ionicProfile = struct.sol(end, :, 3);
acc_zone = x >= p.tp - 2e-7 & x <= p.tp + 10e-7;
x_acc = x(acc_zone) * 1e7;
ionicProfile_acc = ionicProfile(acc_zone);
depl_zone = x >= p.tp + p.ti - 10e-7 & x <= p.tp + p.ti + 2e-7;
x_depl = x(depl_zone) * 1e7;
ionicProfile_depl = ionicProfile(depl_zone);
data_ionic_acc = [x_acc', ionicProfile_acc'];
data_ionic_depl = [x_depl', ionicProfile_depl'];

% copied from pinana
V = struct.sol(:,:,4);     % electric potential
% Calculate energy levels and chemical potential         
V = V - p.EA;                                % Electric potential
Ecb = p.EA-V-p.EA;  % Conduction band potential
n = struct.sol(:,:,1);
P = struct.sol(:,:,2);
Efn = real(-V+p.Ei+(p.kB*p.T/p.q)*log(n/p.ni)); % Electron quasi-Fermi level
Efp = real(-V+p.Ei-(p.kB*p.T/p.q)*log(P/p.ni)); % Hole quasi-Fermi level

Efn = Efn(end, :);
Efp = Efp(end, :);
V = struct.sol(end,:,4);     % electric potential
% Calculate energy levels and chemical potential         
V = V - p.EA;                                % Electric potential
Ecb = p.EA-V-p.EA;                             % Conduction band potential
Evb = p.IP-V-p.EA;
data_band = [xnm', Ecb', Efn', Efp', Evb'];

%% join fields

toBeSavedIonic_acc = [string(headerIonic); string(unitsIonic); string(data_ionic_acc)];

toBeSavedIonic_depl = [string(headerIonic); string(unitsIonic); string(data_ionic_depl)];

toBeSavedBand = [string(headerBand); string(unitsBand); string(data_band)];

%% save csv

fid_ionic_acc = fopen([prefix '-ss-ionic_acc.txt'], 'wt+');
fid_ionic_depl = fopen([prefix '-ss-ionic_depl.txt'], 'wt+');
fid_band = fopen([prefix '-ss-band_diagram.txt'], 'wt+');


for i = 1:size(toBeSavedIonic_acc, 1)
    fprintf(fid_ionic_acc, '%s\t', toBeSavedIonic_acc(i, :));
    fprintf(fid_ionic_acc, '\n');
end

for i = 1:size(toBeSavedIonic_depl, 1)
    fprintf(fid_ionic_depl, '%s\t', toBeSavedIonic_depl(i, :));
    fprintf(fid_ionic_depl, '\n');
end

for i = 1:size(toBeSavedBand, 1)
    fprintf(fid_band, '%s\t', toBeSavedBand(i, :));
    fprintf(fid_band, '\n');
end

fclose(fid_ionic_acc);
fclose(fid_ionic_depl);
fclose(fid_band);
