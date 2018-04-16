function fromISwaveEAStructToTxt(struct, prefix)
% save the main data from a oscilalting solution created by ISwave_EA_single_exec to
% txt files, ideally easy to import with Origin (from OriginLab)

p = struct.p;

%% create header

headerBand = ["x", "CB_Vmax", "Efn_Vmax", "Efp_Vmax", "CB_V0decr", "Efn_V0decr", "Efp_V0decr", "CB_Vmin", "Efn_Vmin", "Efp_Vmin", "CB_V0incr", "Efn_V0incr", "Efp_V0incr"];
headerCurrent = ["t", "Vapp", "total", "ionic", "nonionic", "accumulating"];

%% get measure units

unitsBand = ["nm", repelem("eV", length(headerBand)-1)];
unitsCurrent = ["s", "V", repelem("mA/cm\+(2)", length(headerCurrent)-2)];

%% get data - band diagram

x = struct.x;
xnm = x * 1e7;

t_start = struct.p.tmax - 1.5*pi/struct.p.Vapp_params(end); % last oscillating period
t_array = linspace(t_start, struct.p.tmax, 4);
% pre-allocate
t_index = zeros(length(t_array), 1);
for i = 1:length(t_array)
    [~, t_index(i)] = min(abs(struct.t - t_array(i))); % Finds the index for the time value given
end

[~, ~, ~, Efn, Efp, ~] = pinAna(struct);

V = struct.sol(:,:,4);     % electric potential
% Calculate energy levels and chemical potential         
V = V - p.EA;                                % Electric potential
Ecb = p.EA-V-p.EA;  % Conduction band potential

Ecb_Vmax = Ecb(t_index(1), :);
Efn_Vmax = Efn(t_index(1), :);
Efp_Vmax = Efp(t_index(1), :);

Ecb_V0decr = Ecb(t_index(2), :);
Efn_V0decr = Efn(t_index(2), :);
Efp_V0decr = Efp(t_index(2), :);

Ecb_Vmin = Ecb(t_index(3), :);
Efn_Vmin = Efn(t_index(3), :);
Efp_Vmin = Efp(t_index(3), :);

Ecb_V0incr = Ecb(t_index(4), :);
Efn_V0incr = Efn(t_index(4), :);
Efp_V0incr = Efp(t_index(4), :);

data_band = [xnm', Ecb_Vmax', Efn_Vmax', Efp_Vmax', Ecb_V0decr', Efn_V0decr', Efp_V0decr', Ecb_Vmin', Efn_Vmin', Efp_Vmin', Ecb_V0incr', Efn_V0incr', Efp_V0incr'];

%% get data - currents

time = struct.t; % row

Vapp = p.Vapp_func(p.Vapp_params, struct.t); % row

Jn = struct.Jn; % in mA, column

% block taken from ISwave_single_analysis
% Intrinsic points logical array
itype_points= (struct.x >= p.tp & struct.x <= p.tp + p.ti);
% subtract background ion concentration, for having less noise in trapz
i_matrix = struct.sol(:, :, 3) - p.NI;
% calculate electric field due to ions
Efield_i = p.e * cumtrapz(struct.x, i_matrix, 2) / p.eppi;
% an average would be enough if the spatial mesh was homogeneous in the
% intrinsic, indeed I have to use trapz for considering the spatial mesh
Efield_i_mean = trapz(struct.x(itype_points), Efield_i(:, itype_points), 2) / p.ti;
% calculate displacement current due to ions
Ji_disp = -p.eppi * gradient(Efield_i_mean, struct.t); % in Amperes, column

J_noionic = Jn - Ji_disp * 1000; % in mA, column

% try to calculate the accumulated charge
[dQ_t, ~, ~, ~, ~] = ISwave_subtracting_analysis(struct);

J_accumulating = -dQ_t * 1000; % mA, row

data_current = [time', Vapp', Jn, Ji_disp, J_noionic, J_accumulating'];

%% join fields

toBeSavedBand = [string(headerBand); string(unitsBand); string(data_band)];
toBeSavedCurrent = [string(headerCurrent); string(unitsCurrent); string(data_current)];

%% save csv

fid_band = fopen([prefix '-ac_band_diagram.txt'], 'wt+');
fid_current = fopen([prefix '-ac_current.txt'], 'wt+');

for i = 1:size(toBeSavedBand, 1)
    fprintf(fid_band, '%s\t', toBeSavedBand(i, :));
    fprintf(fid_band, '\n');
end

for i = 1:size(toBeSavedCurrent, 1)
    fprintf(fid_current, '%s\t', toBeSavedCurrent(i, :));
    fprintf(fid_current, '\n');
end

fclose(fid_band);
fclose(fid_current);
