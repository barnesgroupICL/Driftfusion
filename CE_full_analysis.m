function CE_full_analysis(CE_struct)
%CE_FULL_ANALYSIS - Plot Charge Extraction (CE) in a range of background light intensity versus open circuit voltage
%
% Syntax:  CE_full_analysis(CE_struct)
%
% Inputs:
%   CE_STRUCT - a struct containing the most important results of the CE simulation
%
% Example:
%   CE_full_analysis(CE_struct)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also CE_full_exec.

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
set(0, 'defaultLineLineWidth', 2);

F = @(coef, xdata) coef(1) * xdata + coef(2) * (exp(xdata * coef(3)) - 1);

figure('Name', ['CE at light intensity - comparing with reference - ' CE_struct.sol_name], 'NumberTitle', 'off')
    hold off
    plot(CE_struct.Voc, CE_struct.extracting_charges, 'r+'); % extracted by CE
    hold on
    plot(CE_struct.Voc, CE_struct.early_extracting_charges, 'rx'); % just considering first time part of CE
    plot(CE_struct.Voc, CE_struct.subtracting_charges, 'bs'); % calculated delta electron density
    plot(CE_struct.Voc, CE_struct.subtracting_charges_contacts, 'b+'); % calculated from electron density in the contacts
    plot(CE_struct.Voc, CE_struct.subtracting_charges_intrinsic, 'bx'); % calculated from electron density in the intrinsic
%    plot(CE_struct.Voc, CE_struct.ions_displaced, 'bd'); % calculated from ionic delta concentration
    
    legend_array = ["Extracted Charge", "Early Extr Charge", "Subtracting Charge Profiles", "Subtracting Charge at Contacts", "Subtracting Charge in Intrinsic"];
       
    % at the end for not including in the legend
    plot(CE_struct.Voc, F([CE_struct.extracting_geom_cap, CE_struct.extracting_chem_cap_amp, CE_struct.extracting_chem_cap_gamma], CE_struct.Voc), 'r');
    plot(CE_struct.Voc, F([CE_struct.early_geom_cap, CE_struct.early_chem_cap_amp, CE_struct.early_chem_cap_gamma], CE_struct.Voc), 'r');
    plot(CE_struct.Voc, F([CE_struct.subtracting_intrinsic_geom_cap, CE_struct.subtracting_intrinsic_chem_cap_amp, CE_struct.subtracting_intrinsic_chem_cap_gamma], CE_struct.Voc), 'b');

    xlabel('Voc due to light bias [V]');
    ylabel('Extracted Charge [C/cm^2]');
    legend(legend_array);
    
    
figure('Name', ['CE at light intensity removing geometric - comparing with reference - ' CE_struct.sol_name], 'NumberTitle', 'off')
    hold off
    plot(CE_struct.Voc, CE_struct.extracting_charges - F([CE_struct.extracting_geom_cap, 0, 0], CE_struct.Voc), 'r+'); % extracted by CE
    hold on
    plot(CE_struct.Voc, CE_struct.early_extracting_charges - F([CE_struct.early_geom_cap, 0, 0], CE_struct.Voc), 'rx'); % just considering first time part of CE
    plot(CE_struct.Voc, CE_struct.subtracting_charges_intrinsic, 'bx'); % calculated from electron density in the intrinsic
%    plot(CE_struct.Voc, CE_struct.ions_displaced, 'bd'); % calculated from ionic delta concentration
    
    legend_array2 = ["Extracted Charge nogeom", "Early Extr Charge nogeom", "Subtracting Charge in Intrinsic"];
       
    % at the end for not including in the legend
    plot(CE_struct.Voc, F([0, CE_struct.extracting_chem_cap_amp, CE_struct.extracting_chem_cap_gamma], CE_struct.Voc), 'r');
    plot(CE_struct.Voc, F([0, CE_struct.early_chem_cap_amp, CE_struct.early_chem_cap_gamma], CE_struct.Voc), 'r');
    plot(CE_struct.Voc, F([0, CE_struct.subtracting_intrinsic_chem_cap_amp, CE_struct.subtracting_intrinsic_chem_cap_gamma], CE_struct.Voc), 'b');

    xlabel('Voc due to light bias [V]');
    ylabel('Extracted Charge [C/cm^2]');
    legend(legend_array2);
    
    
%------------- END OF CODE --------------