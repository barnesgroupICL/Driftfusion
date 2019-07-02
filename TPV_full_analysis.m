function TPV_full_analysis(TPV_struct)
%TPV_FULL_ANALYSIS - Plot Transient PhotoVoltage (TPV) in a range of background light intensity versus open circuit voltage
% Both mono and bi exponential fitting results are plotted versus
% light-bias
%
% Syntax:  TPV_full_analysis(TPV_struct)
%
% Inputs:
%   TPV_STRUCT - a struct containing the most important results of the TPV simulation
%
% Example:
%   TPV_full_analysis(TPV_struct)
%     do plot
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also TPVvariab_full_exec, TPVconst_full_exec.

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

figure('Name', ['Mono-exponential fit TPV ' TPV_struct.sol_name], 'NumberTitle', 'off')
    yyaxis left
    hold off
    semilogy(TPV_struct.Voc, TPV_struct.monoexp_tau, 'bo-');
    xlabel('Voc [V]');
    ylabel('Tau [s]');
    yyaxis right
    hold off
    plot(TPV_struct.Voc, TPV_struct.monoexp_lin, 'ro');
    ylabel('y0 [V]');
    legend('Tau', 'y0');

figure('Name', ['Bi-exponential fit TPV ' TPV_struct.sol_name], 'NumberTitle', 'off')
    yyaxis left
    hold off
    semilogy(TPV_struct.Voc, TPV_struct.biexp_tau_fast, 'b.-');
    hold on
    semilogy(TPV_struct.Voc, TPV_struct.biexp_tau_slow, 'bo-');
    xlabel('Voc [V]');
    ylabel('Tau [s]');
    hold off
    yyaxis right
    hold off
    plot(TPV_struct.Voc, TPV_struct.biexp_lin_fast, 'r.');
    hold on
    plot(TPV_struct.Voc, TPV_struct.biexp_lin_slow, 'ro');
    ylabel('y0 [V]');
    legend('Tau fast', 'Tau slow', 'y0 fast', 'y0 slow');
    hold off
    
%------------- END OF CODE --------------