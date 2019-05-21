function plot_charges_single(struct, images_prefix)
%PLOT_CHARGES_SINGLE - Plot free charges and ions densities, and energy levels over device thickness and over time.
% This is inspired by Phil Calado's pinExtract2.m but rewritten
%
% Syntax:  plot_charges_single(struct, save_name)
%
% Inputs:
%   STRUCT - a symmetric or asymmetric struct of a solution evolving over
%     time
%   IMAGES_PREFIX - optional char array, if provided saves the images in the
%     current directory with the name starting with images_prefix
%
% Example:
%   plot_charges_single(TPV_single_exec(ssol_i_light, 1e-3, 5))
%     do plot of TPV decay, do not save images
%   plot_charges_single(CE_single_exec(asymmetricize(ssol_i_light), 1e-3, 1), 'CEexp-')
%     do plot of CE experiment and save images
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

% increase graphics font size
set(0, 'defaultAxesFontSize', 24);
% set image dimension
set(0, 'defaultfigureposition', [0, 0, 1000, 750]);
% set line thickness
set(0, 'defaultLineLineWidth', 1);

% criteria for recognizing impedance spectroscopy
is_t_linear = struct.p.tmesh_type == 1;

%% define times and positions interesting for plotting

% Array containing times [s] for which you want to extract the data
t_array = [1e-10]; % in seconds, gets expanded to longer times further in the code with a logspace function

% shortcuts
tp = 1e7*struct.p.tp;
ti = 1e7*struct.p.ti;

if struct.p.mui
    % Array containing position [nm] for which you want to extract the data
%    x_array = [0, tp, tp, tp + ti/2, tp + ti/2, tp + ti/2, tp + ti/2, tp + ti, tp + ti, 0]; % last and first values gets obtained by the maximum conc position of p and n
    x_array = [tp-2, tp, tp+2, tp, tp+ti, tp+ti-2, tp+ti, tp+ti+2]; % last and first values gets obtained by the maximum conc position of p and n

    % select a specific kind of charge to be plotted at each spatial point
%    x_ion_boolean =         logical([0, 1, 0, 1, 0, 0, 0, 1, 0, 0]);
%    x_electron_boolean =    logical([0, 0, 1, 0, 1, 0, 0, 0, 0, 1]);
%    x_hole_boolean =        logical([1, 0, 0, 0, 0, 1, 0, 0, 1, 0]);
   x_ion_boolean =         logical([0, 0, 0, 1, 1, 0, 0, 0]);
   x_electron_boolean =    logical([1, 1, 1, 0, 0, 0, 0, 0]);
   x_hole_boolean =        logical([0, 0, 0, 0, 0, 1, 1, 1]);

    % curves are not normalized individually, but with
    % respect to the maximum absolute value of taken from group of curves
    % if groups have to space over different kind of plotted charges more
    % code is needed, up to now goups works just for ions with ions,
    % electrons with holes, raw charge with raw charge
%    x_normalize_groups = ["np_contact", "a_interf", "np_contact", "a_intr", "np_intr", "np_intr", "net_intr", "a_interf", "np_contact", "np_contact"];
    x_normalize_groups = ["np_contact", "np_contact", "np_contact", "a_interf", "a_interf", "np_contact", "np_contact", "np_contact"];

else
    x_array = [0, 400, 400, 0]; % last and first values gets obtained by the maximum conc position of p and n
    x_ion_boolean = false(length(x_array), 1);
    x_electron_boolean =    logical([0, 1, 0, 1]);
    x_hole_boolean =        logical([1, 0, 1, 0]);
    x_normalize_groups = [1, 2, 2, 1];
end

t_points = length(struct.sol(:, 1, 1)); % this works even if the computation was not completed
x_points = length(struct.sol(1, :, 1)); % struct.p.xpoints is 0 mysteriously, so can't be used

%% define largest time and time array

if ~is_t_linear % in Impedance Spectroscopy simulation tmesh_type is 1 and the time mesh is linear
    if isfield(struct, 'Voc') % if solution is symmetric
        [deltaV, t_max_Voc_index] = max(struct.Voc);
        [~, t_index_tenth_Voc] = min(abs(struct.Voc(t_max_Voc_index:end) - 0.1 * deltaV));
        t_tenth_Voc = struct.t(min(t_index_tenth_Voc + t_max_Voc_index, length(struct.t)));
        t_quasiend = t_tenth_Voc;
    else
        [~, t_index_tenth_tend] = min(abs(struct.t - 0.1 * struct.t(t_points)));
        t_quasiend = struct.t(t_index_tenth_tend);
    end

    t_quasiend = max(t_quasiend, 1e-2 * struct.t(t_points));
    t_array = [t_array, logspace(log10(max(t_array)), log10(t_quasiend) + 0.1, 7)]; % duplicates gets eliminated further in the code
else
    % struct.p.Vapp_params[end] is 2pi*frequency of ISwave
    t_start = struct.p.tmax - 1.5*pi/struct.p.Vapp_params(end); % last oscillating period
    t_array = linspace(t_start, struct.p.tmax, 4);
    t_quasiend = t_points(end); % just needed for graphing
end

t_array = sort(t_array);

%% data used by both plots over space and over time

xnm = struct.x*1e7; % position in [nm]

% extract particles matrix
n = struct.sol(:, :, 1);
p = struct.sol(:, :, 2);
a = struct.sol(:, :, 3);

% Electric potential, EA usually is zero
V = struct.sol(:, :, 4) - struct.p.EA;
% Conduction band potential
Ecb = -V;
% Valence band potential
Evb = struct.p.IP - V - struct.p.EA;

% p-type binary matrix
pBM = ones(t_points, x_points)*diag(struct.x <= struct.p.tp);
% Intrinsic binary matrix
iBM = ones(t_points, x_points)*diag(struct.x > struct.p.tp & struct.x < struct.p.tp + struct.p.ti);
% n-type binary matrix
nBM = ones(t_points, x_points)*diag(struct.x >= struct.p.tp + struct.p.ti & struct.x <= struct.p.tp + struct.p.ti + struct.p.tn);

nstat = (- struct.p.NA - struct.p.NI)*pBM + (-struct.p.NI)*iBM + (struct.p.ND-struct.p.NI)*nBM;

if struct.p.OC
    nstat = nstat(:, 1:floor(x_points / 2) + 1); % take just the first spatial half
    % symmetricize in the same way as symmetricize.m
    nstat2 = fliplr(nstat(:,:,:));
    nstat2 = nstat2(:, 2:end, :); % Delete middle point
    nstat = [nstat, nstat2]; % concatenate
end
rhoc = (-n + p + a + nstat); % Net charge density calculated from adding individual charge densities

% reconstruct Efn and Efp as in pindrift
Ei = struct.p.Ei; % Intrinsic Fermi Energy
kB = struct.p.kB;
T = struct.p.T;
q = struct.p.q;
ni = struct.p.ni; % Intrinsic carrier density
Efn = real(-V+Ei+(kB*T/q)*log(n/ni)); % Electron quasi-Fermi level
Efp = real(-V+Ei-(kB*T/q)*log(p/ni)); % Hole quasi-Fermi level

%% build t_array and data for plotting in specific times

% use this for not showing the full device in case of symmetrical device
if struct.p.OC
    xnm_end = xnm(ceil(end/2));
else
    xnm_end = xnm(end);
end

% pre-allocate
t_index = zeros(length(t_array), 1);
t_legend_text = string(t_index);

for i = 1:length(t_array)
    [~, t_index(i)] = min(abs(struct.t - t_array(i))); % Finds the index for the time value given
    t_legend_text(i) = strcat(num2str(t_array(i)), ' s');
end

% remove time points out of the maximum time resulting in duplicate lines
[t_index, t_unique] = unique(t_index);
t_array = t_array(t_unique);

% in case of Impedance Spectroscopy, throw away previous legend a set a
% more meaningful one
if ~is_t_linear
    t_legend_text = t_legend_text(t_unique);
else
    t_legend_text = ["max V AC peak", "zero V AC decreasing", "min V AC peak", "zero V AC increasing"];
end

%% common settings for markers in plots

markers = {'o', 'x'}; % for plotting markers on top of lines
mark_spac = 60; % points between one marker and the next one
mark_star = round(mark_spac/length(t_array)); % parameter used for position of first marker for each line

%% find the position of the peak of free charges in the contacts

% these picks also peaks in the "wrong" direction, so even if the charge
% delta is not an increase of mayor carrier for that region
ptype_peak_overtime = max(abs(struct.sol(t_index, logical(pBM(1,:)), 2) - struct.sol(1, logical(pBM(1,:)), 2))); % results in a spacial array, the time dimension gets eliminated byt the max function
[ptype_y_limit, peak_ptype] = max(abs(ptype_peak_overtime)); % just the first half of a symmetrical solution is present
if ~x_array(1)
    x_array(1) = 1e7 * struct.x(peak_ptype);
end

ntype_peak_overtime = max(abs(struct.sol(t_index, logical(nBM(1,:)), 1) - struct.sol(1, logical(nBM(1,:)), 1))); % results in a spacial array, the time dimension gets eliminated byt the max function
[~, peak_ntype] = max(abs(ntype_peak_overtime)); % just the first half of a symmetrical solution is present
if ~x_array(end)
    x_array(end) = 1e7 * struct.x(nnz(pBM(1,:)) + nnz(iBM(1,:)) + peak_ntype);
end

%% build x_array and data for plotting in specific positions

% pre-allocate
x_index = zeros(length(x_array), 1);
rhoc_selected_deltas = zeros(length(x_array), t_points);
a_selected_deltas = rhoc_selected_deltas;
n_selected_deltas = rhoc_selected_deltas;
p_selected_deltas = rhoc_selected_deltas;
normalize_rhoc_array = zeros(length(x_array), 1);
normalize_a_array = normalize_rhoc_array;
normalize_np_array = normalize_rhoc_array;

for i = 1:length(x_array)
    % Finds the index for the position value given
    [~, x_index(i)] = min(abs(xnm - x_array(i)));
    % get delta charge profiles at specific positions
    rhoc_selected_deltas(i, :) = transpose(rhoc(:, x_index(i)) - rhoc(1, x_index(i)));
    % get delta ionic charge profiles at specific positions
    if struct.p.mui
        a_selected_deltas(i, :) = transpose(a(:, x_index(i)) - a(1, x_index(i)));
    end
    % get delta electronic charge profiles at specific positions
    n_selected_deltas(i, :) = transpose(n(:, x_index(i)) - n(1, x_index(i)));
    % get delta holes charge profiles at specific positions
    p_selected_deltas(i, :) = transpose(p(:, x_index(i)) - p(1, x_index(i)));
end

% loop again for normalization
for i = 1:length(x_array)
    % find indexes of a group of curves to be normalized together
    normalize_array = x_normalize_groups == x_normalize_groups(i);
    normalize_rhoc_array(i) = max(max(abs(rhoc_selected_deltas(normalize_array, :))));
    normalize_a_array(i) = max(max(abs(a_selected_deltas(normalize_array, :))));
    normalize_np_array(i) = max(max(max(abs(n_selected_deltas(normalize_array, :)))), max(max(abs(p_selected_deltas(normalize_array, :)))));
end

% divide each row in each matrix by the correct normalization factor
rhoc_selected_deltas = rhoc_selected_deltas ./ normalize_rhoc_array; % normalize
a_selected_deltas = a_selected_deltas ./ normalize_a_array; % normalize
n_selected_deltas = n_selected_deltas ./ normalize_np_array; % normalize
p_selected_deltas = p_selected_deltas ./ normalize_np_array; % normalize

%% obtain total charge, taken from CE_ISstep_subtracting_analysis

n_delta_profile = struct.sol(:, :, 1) - struct.sol(1, :, 1);

% sum of the electrons as they could be extracted at the
% respective electrode completely
n_delta = trapz(struct.x, n_delta_profile, 2);
subtracting_charges = -struct.p.e * n_delta; % convert from number of electrons to charge
subtracting_charges_norm = subtracting_charges / range(subtracting_charges);

%% plots at specific times over all positions

% create a color array with one color more than necessary
if length(t_array) == 5 % too many blue lines in this case
    jet_matrix = jet(length(t_array) + 2);
    jet_matrix = jet_matrix([1,3,4,5,6,7], :);
else
    jet_matrix = jet(length(t_array) + 1);
end
% find the yellow (which in RGB code is 1,1,0) and remove it from the
% colors list
jet_yellow_logical = ismember(jet_matrix, [1, 1, 0], 'rows');
jet_no_yellow = jet_matrix(~jet_yellow_logical, :);
jet_no_yellow_flip = flipud(jet_no_yellow);
t_colors = colormap(jet_no_yellow_flip);

figure('Name', 'Net Charge Density at specific times [s]', 'NumberTitle', 'off')
    hold on
    yyaxis right
    ylabel('V [V]') % right y-axis
    for i = 1:length(t_index)
        plot(xnm, Ecb(t_index(i), :), 'Color', t_colors(i, :), 'Marker', 'none', 'LineStyle', '--');
    end
    
    yyaxis left
    for i = 1:length(t_index)
        plot(xnm, rhoc(t_index(i), :), 'Color', t_colors(i, :), 'Marker', 'none', 'LineStyle', '-');
    end
    
    xlabel('Position [nm]');
    xlim([0, xnm_end]); % in case of symmetrical device this can limit the graph to the first half
    ylabel({'Net Charge';'Density [cm^{-3}]'}); % left y-axis
    legend(t_legend_text);
    % just for verifying that the signs were correct, derive charge from
    % potential second derivative
%     charge = 1e21 * gradient(gradient(Ecb(t_index(2), :))./gradient(xnm))./gradient(xnm);
%     plot(xnm, charge, 'LineWidth', 3, 'LineStyle', '-');
    hold off

    if exist('images_prefix', 'var')
        fig = gcf;
        saveas(fig, strcat(string(images_prefix), '-charge_vs_position.png'));
        saveas(fig, strcat(string(images_prefix), '-charge_vs_position.fig'));
    end

he = zeros(length(t_index), 1);
figure('Name', 'Delta Net Charge Density and Energy Levels at specific times [s]', 'NumberTitle', 'off')
    hold on
    
    yyaxis left
    for i = 1:length(t_index)
        plot(xnm, rhoc(t_index(i), :) - rhoc(1, :), 'Color', t_colors(i, :), 'Marker', 'none', 'LineStyle', '--', 'LineWidth', 2);
    end
    current_ylim = get(gca,'ylim');
    plot(1e7*[struct.p.tp, struct.p.tp], current_ylim, 'k-'); % draw vertical line for p-type - intrinsic boundary
    plot(1e7*[struct.p.tp + struct.p.ti, struct.p.tp + struct.p.ti], current_ylim, 'k-'); % draw vertical line for intrinsic - n-type boundary
    ylabel({'Delta Charge';'Density [cm^{-3}]'}); % left y-axis

    yyaxis right
    ylabel({'Conduction band';'variation [eV]'}) % right y-axis
    for i = 1:length(t_index)
        he(i) = plot(xnm, Ecb(t_index(i), :) - Ecb(1, :), 'Color', t_colors(i, :), 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 2);
    end
    legend(he, t_legend_text);
    
    xlabel('Position [nm]');
    xlim([0, xnm_end]); % in case of symmetrical device this can limit the graph to the first half
    hold off
    
    if exist('images_prefix', 'var')
        fig = gcf;
        saveas(fig, strcat(string(images_prefix), '-delta_charge_vs_position.png'));
        saveas(fig, strcat(string(images_prefix), '-delta_charge_vs_position.fig'));
        ax = gca;
        ax.XLim = [x_array(1)-40, x_array(1)+40];
        ax.YLim = [-ptype_y_limit, ptype_y_limit];
        saveas(fig, strcat(string(images_prefix), '-delta_charge_vs_position-zoom_ptype.png'));
        saveas(fig, strcat(string(images_prefix), '-delta_charge_vs_position-zoom_ptype.fig'));
        if struct.p.mui
            ax.XLim = [x_array(2)-5, x_array(2)+10];
            peak2 = max(abs(rhoc(:, x_index(2)) - rhoc(1, x_index(2))));
            ax.YLim = [-peak2, peak2];
            saveas(fig, strcat(string(images_prefix), '-delta_charge_vs_position-zoom_intr.png'));
            saveas(fig, strcat(string(images_prefix), '-delta_charge_vs_position-zoom_intr.fig'));
        end
    end
    
figure('Name', 'Delta Electron (dashed), Hole (dot dashed) and Ion (solid) charge at specific times [s]', 'NumberTitle', 'off')
    hold on
    %yyaxis right
    %ylabel('V [V]') % right y-axis
    %for i = 1:length(t_index)
        %plot(xnm, Ecb(t_index(i), :) - Ecb(1, :), 'Color', t_colors(i, :), 'Marker', 'none', 'LineStyle', '--', 'LineWidth', 0.5);
        %plot(xnm(i*mark_star:mark_spac:end), Ecb(t_index(i), i*mark_star:mark_spac:end) - Ecb (1, i*mark_star:mark_spac:end), 'Color', t_colors(i, :), 'MarkerFaceColor', t_colors(i, :), 'LineStyle', 'none', 'Marker', markers{mod(i, numel(markers))+1}, 'MarkerSize', 5);
    %end
    
    %yyaxis left
    Hp = zeros(length(t_index), 1); % pre allocate plot handle for holes
    for i = 1:length(t_index)
        plot(xnm, - n(t_index(i), :) + n(1, :), 'Color', t_colors(i, :), 'Marker', 'none', 'LineStyle', '--', 'LineWidth', 2);
        %plot(xnm(i*mark_star:mark_spac:end), - n(t_index(i), i*mark_star:mark_spac:end) + n(1, i*mark_star:mark_spac:end), 'Color', t_colors(i, :), 'MarkerFaceColor', t_colors(i, :), 'LineStyle', 'none', 'Marker', markers{mod(i, numel(markers))+1}, 'MarkerSize', 5);

        Hp(i) = plot(xnm, p(t_index(i), :) - p(1, :), 'Color', t_colors(i, :), 'Marker', 'none', 'LineStyle', '-.', 'LineWidth', 2);
        %plot(xnm(i*mark_star:mark_spac:end), p(t_index(i), i*mark_star:mark_spac:end) - p(1, i*mark_star:mark_spac:end), 'Color', t_colors(i, :), 'MarkerFaceColor', t_colors(i, :), 'LineStyle', 'none', 'Marker', markers{mod(i, numel(markers))+1}, 'MarkerSize', 5);

        if struct.p.mui % in this case overwrite handles
            Hp(i) = plot(xnm, a(t_index(i), :) - a(1, :), 'Color', t_colors(i, :), 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 2);
            %plot(xnm(i*mark_star:mark_spac:end), a(t_index(i), i*mark_star:mark_spac:end) - a(1, i*mark_star:mark_spac:end), 'Color', t_colors(i, :), 'MarkerFaceColor', t_colors(i, :), 'LineStyle', 'none', 'Marker', markers{mod(i, numel(markers))+1}, 'MarkerSize', 5);    
        end
    end
    current_ylim = get(gca,'ylim');
    plot(1e7*[struct.p.tp, struct.p.tp], current_ylim, 'k-', 'LineWidth', 0.5); % draw vertical line for p-type - intrinsic boundary
    plot(1e7*[struct.p.tp + struct.p.ti, struct.p.tp + struct.p.ti], current_ylim, 'k-', 'LineWidth', 0.5); % draw vertical line for intrinsic - n-type boundary
    
    xlabel('Position [nm]');
    xlim([0, xnm_end]); % in case of symmetrical device this can limit the graph to the first half
    ylabel({'Delta Charge'; 'Density [cm^{-3}]'}); % left y-axis
    legend(Hp, t_legend_text); % create the legend for holes, they have solid lines which are better for legend
    legend boxoff
    hold off
    
    if exist('images_prefix', 'var')
        fig = gcf;
        saveas(fig, strcat(string(images_prefix), '-delta_n_p_a_vs_position.png'));
        saveas(fig, strcat(string(images_prefix), '-delta_n_p_a_vs_position.fig'));
        ax = gca;
        ax.XLim = [x_array(1)-40, x_array(1)+40];
        ax.YLim = [-ptype_y_limit, ptype_y_limit];
        saveas(fig, strcat(string(images_prefix), '-delta_n_p_a_vs_position-zoom_ptype.png'));
        saveas(fig, strcat(string(images_prefix), '-delta_n_p_a_vs_position-zoom_ptype.fig'));
        if struct.p.mui
            ax.XLim = [x_array(2)-5, x_array(2)+10];
            peak2 = max(abs(rhoc(:, x_index(2)) - rhoc(1, x_index(2))));
            ax.YLim = [-peak2, peak2];
            saveas(fig, strcat(string(images_prefix), '-delta_n_p_a_vs_position-zoom_intr.png'));
            saveas(fig, strcat(string(images_prefix), '-delta_n_p_a_vs_position-zoom_intr.fig'));
        end
    end
    
h = zeros(length(t_index), 1); 
figure('Name', 'Band Diagram and Quasi Fermi Levels', 'NumberTitle', 'off');
    hold on
    for i = 1:length(t_index)
        plot(xnm, Efn(t_index(i), :), 'Color', t_colors(i, :), 'LineStyle', '--'); % Electrons Quasi Fermi
        plot(xnm, Efp(t_index(i), :), 'Color', t_colors(i, :), 'LineStyle', '--'); % Holes Quasi Fermi
        plot(xnm, Ecb(t_index(i), :), 'Color', t_colors(i, :), 'LineStyle', '-'); % Conduction Band
        plot(xnm, Evb(t_index(i), :), 'Color', t_colors(i, :), 'LineStyle', '-'); % Valence Band
    end
    % separated for cycles for having markers drawn always on top of lines
    % stack, useful for overlapping lines covering the markers
    for i = 1:length(t_index)
        plot(xnm(i*mark_star:mark_spac:end), Efn(t_index(i), i*mark_star:mark_spac:end), 'Color', t_colors(i, :), 'MarkerFaceColor', t_colors(i, :), 'LineStyle', 'none', 'Marker', markers{mod(i, numel(markers))+1}, 'MarkerSize', 5); % Electrons Quasi Fermi
        plot(xnm(i*mark_star:mark_spac:end), Efp(t_index(i), i*mark_star:mark_spac:end), 'Color', t_colors(i, :), 'MarkerFaceColor', t_colors(i, :), 'LineStyle', 'none', 'Marker', markers{mod(i, numel(markers))+1}, 'MarkerSize', 5); % Holes Quasi Fermi
        h(i) = plot(xnm(i*mark_star:mark_spac:end), Ecb(t_index(i), i*mark_star:mark_spac:end), 'Color', t_colors(i, :), 'MarkerFaceColor', t_colors(i, :), 'LineStyle', 'none', 'Marker', markers{mod(i, numel(markers))+1}, 'MarkerSize', 5); % Conduction Band
        plot(xnm(i*mark_star:mark_spac:end), Evb(t_index(i), i*mark_star:mark_spac:end), 'Color', t_colors(i, :), 'MarkerFaceColor', t_colors(i, :), 'LineStyle', 'none', 'Marker', markers{mod(i, numel(markers))+1}, 'MarkerSize', 5); % Valence Band
    end
    range_ylim = range(get(gca, 'ylim'));
    text(0, max(max(Ecb(:, 1:end*0.3))) + range_ylim * 0.05, '   Conduction band', 'FontSize', 20)
    text(0, max(max(Efn(:, 1:end*0.3))) + range_ylim * 0.05, '   e^- quasi Fermi', 'FontSize', 20)
    text(0, max(max(Efp(:, 1:end*0.3))) + range_ylim * 0.05, '   h^+ quasi Fermi', 'FontSize', 20)
    text(0, min(min(Evb(:, 1:end*0.3))) - range_ylim * 0.05, '   Valence band', 'FontSize', 20)

    ylabel('Energy [eV]');
    xlim([0, xnm_end]);
    xlabel('Position [nm]');
    legend(h, t_legend_text);
    legend boxoff; % removes background and outline from legend

hd = zeros(length(t_index), 1); 
figure('Name', 'Electrons and holes charge densities', 'NumberTitle', 'off');
    hold on
    for i = 1:length(t_index)
        hd(i) = plot(xnm, n(t_index(i), :), 'Color', t_colors(i, :));
        plot(xnm, p(t_index(i), :), '--', 'Color', t_colors(i, :));
    end
    ax = gca;
    ax.YScale = 'log';
    ylabel('{\itn, p} [cm^{-3}]')
    xlim([0, xnm_end]);
    ylim([1e0, 1e20]);
    legend(hd, t_legend_text);

%% plots in specific points over all times

% create a color array with one color more than necessary
jet_matrix = jet(length(x_array) + 1);
% find the yellow (which in RGB code is 1,1,0) and remove it from the
% colors list
jet_yellow_logical = ismember(jet_matrix, [1, 1, 0], 'rows');
jet_no_yellow = jet_matrix(~jet_yellow_logical, :);
jet_no_yellow_flip = flipud(jet_no_yellow);
x_colors = colormap(jet_no_yellow_flip);

% pre-allocate
x_legend_text = string(zeros(length(x_array), 1));

figure('Name', 'Normalized Delta Charge in specific points [nm]', 'NumberTitle', 'off')
    hold on
    yyaxis right
    % if OC or linear time mesh, like in impedance spectroscopy simulation
    % where a linear time mesh is used and voltage is applied
    if is_t_linear
        ylabel('Vapp [V]') % right y-axis
        h(1) = plot(struct.t(1:t_points), struct.p.Vapp_func(struct.p.Vapp_params, struct.t(1:t_points)), 'Marker', 'none');     
        x_legend_text(1) = 'Applied Voltage';
    elseif isfield(struct, 'Voc') % exist('struct.Voc', 'var') does not work
        ylabel('V [V]') % right y-axis
        h(1) = plot(struct.t(1:t_points), struct.Voc, 'Marker', 'none');
        x_legend_text(1) = 'Voc';
    else
        ylabel('J [mA/cm2]') % right y-axis
        h(1) = plot(struct.t(1:t_points), struct.Jn, 'Marker', 'none');
        x_legend_text(1) = 'Current';
    end
    ax = gca;
    if ~is_t_linear
        ax.XScale = 'log'; % for putting the scale in log
    end

    yyaxis left
    set(gca, 'YTickLabel', []); % remove numbering from left axis
    for i = 1:length(x_index)
        if x_ion_boolean(i)
            x_legend_text(i+1) = strcat(num2str(x_array(i)), ' nm ionic charge');
            h(i+1) = plot(struct.t(1:t_points), a_selected_deltas(i, :), 'Color', x_colors(i, :), 'Marker', 'o', 'LineWidth', 2, 'LineStyle', '-');
        elseif x_hole_boolean(i)
            x_legend_text(i+1) = strcat(num2str(x_array(i)), ' nm holes charge');
            h(i+1) = plot(struct.t(1:t_points), p_selected_deltas(i, :), 'Color', x_colors(i, :), 'Marker', 'none', 'LineWidth', 2, 'LineStyle', '-');
        elseif x_electron_boolean(i)
            x_legend_text(i+1) = strcat(num2str(x_array(i)), ' nm electronic charge');
            h(i+1) = plot(struct.t(1:t_points), -n_selected_deltas(i, :), 'Color', x_colors(i, :), 'Marker', 'none', 'LineWidth', 2, 'LineStyle', '-');
        else
            x_legend_text(i+1) = strcat(num2str(x_array(i)), ' nm net charge');
            h(i+1) = plot(struct.t(1:t_points), rhoc_selected_deltas(i, :), 'Color', x_colors(i, :), 'Marker', 'none', 'LineWidth', 2, 'LineStyle', '-');
        end
    end
    if is_t_linear && exist('struct.Jn', 'var') % current was not plotted, add it normalized
        J_nobias = struct.Jn - mean([max(struct.Jn), min(struct.Jn)]);
        h(length(h)+1) = plot(struct.t(1:t_points), -2*J_nobias/range(J_nobias), 'b--', 'Marker', 'none');
        x_legend_text = [x_legend_text; 'Current'];
    end
    h(length(h)+1) = plot(struct.t(1:t_points), subtracting_charges_norm, 'k-', 'Marker', 'none');
    x_legend_text = [x_legend_text; 'Total charge'];

    xmin = max(t_array(1) / 10, t_quasiend / 1e11);

    xlim([xmin, struct.t(t_points)]); % ignore first points
    ylim([-1.05, 1.05]); % sane limits for normalized plot
    xlabel('Time [s]');
%    ylabel({'Normalized Delta Charge'}); % left y-axis
    legend(h, x_legend_text);%, 'Interpreter','latex'); % LaTeX can be used for plotting ominus oplus 
    
    hold off

    if exist('images_prefix', 'var')
        fig = gcf;
        saveas(fig, strcat(string(images_prefix), '-delta_charge_vs_time.png'));
        saveas(fig, strcat(string(images_prefix), '-delta_charge_vs_time.fig'));
    end
    
%------------- END OF CODE --------------
