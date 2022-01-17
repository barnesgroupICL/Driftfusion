function IS_EA_struct_exporter(prefix, struct)
%IS_EA_STRUCT_EXPORTER - Exports single impedance simulation data to text file
% Saves the main data from a oscillating solution created by doIS_EA to text files
% for easing the import with Origin (from OriginLab).
%
% Syntax:  IS_EA_struct_exporter(prefix, struct)
%
% Inputs:
%   PREFIX - char array, prefix to be used for the text files names
%   STRUCT - a struct with a solution being perturbed by an
%     oscillating voltage, as generated from IS_EA_single_exec
%
% Example:
%   IS_EA_struct_exporter('asymssol_i_1S_SR_is_100mHz_2mV', asymssol_i_1S_SR_is_100mHz_2mV)
%     save single simulation data to text files
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also doIS_EA, IS_script, EA_script.

%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%------------- BEGIN CODE --------------

p = struct.par;

%% create header

headerBand = ["x", "CB_Vmax", "VB_Vmax", "Efn_Vmax", "Efp_Vmax", "CB_V0decr", "VB_V0decr", "Efn_V0decr", "Efp_V0decr", "CB_Vmin", "VB_Vmin", "Efn_Vmin", "Efp_Vmin", "CB_V0incr", "VB_V0incr", "Efn_V0incr", "Efp_V0incr"];
headerCurrent = ["t", "Vapp", "total", "accumulating", "total_inphase", "total_outofphase", "recombination"];

%% get measure units

unitsBand = ["nm", repelem("eV", length(headerBand)-1)];
unitsCurrent = ["s", "V", repelem("mA/cm\+(2)", length(headerCurrent)-2)];

%% get data - band diagram

x = struct.x;
xnm = x * 1e7;

t_start = struct.par.tmax - 0.75/struct.par.V_fun_arg(3); % last oscillating period
t_array = linspace(t_start, struct.par.tmax, 4);
% pre-allocate
t_index = zeros(length(t_array), 1);
for i = 1:length(t_array)
    [~, t_index(i)] = min(abs(struct.t - t_array(i))); % Finds the index for the time value given
end

[Ecb, Evb, Efn, Efp] = dfana.QFLs(struct);

Ecb_Vmax = Ecb(t_index(1), :);
Evb_Vmax = Evb(t_index(1), :);
Efn_Vmax = Efn(t_index(1), :);
Efp_Vmax = Efp(t_index(1), :);

Ecb_V0decr = Ecb(t_index(2), :);
Evb_V0decr = Evb(t_index(2), :);
Efn_V0decr = Efn(t_index(2), :);
Efp_V0decr = Efp(t_index(2), :);

Ecb_Vmin = Ecb(t_index(3), :);
Evb_Vmin = Evb(t_index(3), :);
Efn_Vmin = Efn(t_index(3), :);
Efp_Vmin = Efp(t_index(3), :);

Ecb_V0incr = Ecb(t_index(4), :);
Evb_V0incr = Evb(t_index(4), :);
Efn_V0incr = Efn(t_index(4), :);
Efp_V0incr = Efp(t_index(4), :);

data_band = [xnm', Ecb_Vmax', Evb_Vmax', Efn_Vmax', Efp_Vmax', Ecb_V0decr', Evb_V0decr', Efn_V0decr', Efp_V0decr', Ecb_Vmin', Evb_Vmin', Efn_Vmin', Efp_Vmin', Ecb_V0incr', Evb_V0incr', Efn_V0incr', Efp_V0incr'];

%% get data - currents

time = struct.t; % row

Vapp_fun = fun_gen(p.V_fun_type);
Vapp = Vapp_fun(p.V_fun_arg, struct.t); % row

J = dfana.calcJ(struct); % in mA, column
Jn = J.tot(:,end);

% try to calculate the accumulated charge
[dQ_t, ~] = IS_ana_subtracting(struct);

J_accumulating = -dQ_t * 1000; % mA, row

coeff = IS_ana_plot(struct, false, true);
cJ = coeff.Jtot;
% in phase electronic current. p.V_fun_arg(3) is the applied voltage
% frequency when p.V_fun_type is sin.
Jn_inphase = Vapp_fun([cJ(1)*cos(cJ(3)), cJ(2)*cos(cJ(3)), p.V_fun_arg(3), 0], struct.t);
% out of phase electronic current
Jn_quadrature = Vapp_fun([cJ(1)*sin(cJ(3)), cJ(2)*sin(cJ(3)), p.V_fun_arg(3), pi/2], struct.t);

% recombination current
r = dfana.calcr(struct); % mA

data_current = [time', Vapp', Jn, J_accumulating', 1000*Jn_inphase', 1000*Jn_quadrature', r.tot(:,end)];

if struct.par.mobseti && any(struct.par.mucat)
    headerCurrentIonic = "ionic";%, "nonionic", "nonionic_inphase", "nonionic_outofphase"];
    headerCurrent = [headerCurrent, headerCurrentIonic];
    unitsCurrent = [unitsCurrent, repelem("mA/cm\+(2)", length(headerCurrentIonic))];
    data_current = [data_current, J.c(:,end) + J.a(:,end)];%, J_noionic, 1000*Jn_noionic_inphase', 1000*Jn_noionic_quadrature'];
end


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

%------------- END OF CODE --------------
