function nid_array = ideality_from_dark_jVpoints(asymstruct_eq, Vend, deltaV)
%IDEALITY_FROM_DARK_JVPOINTS - Obtains voltage dependent ideality factor from the stabilised current at different voltages
% A voltage dependent ideality factor can be obtained from a dark current-voltage sweep
% (or stabilised points in the case of hysteretic devices)
% https://www.pveducation.org/pvcdrom/characterisation/measuring-ideality-factor
%
% Syntax:  nid_array = ideality_from_dark_jVpoints(asymstruct_eq, Vend, deltaV)
%
% Inputs:
%   ASYMSTRUCT_EQ - an asymmetric solution structure as created by PINDRIFT in dark conditions;
%   VEND - the maximum voltage point;
%   DELTAV - the distance between voltage points.
%
% Outputs:
%   NID_ARRAY - an array of ideality factor at different applied voltages.
%
% Example:
%   nid_array = ideality_from_dark_jVpoints(sol_i_eq_SR, 1, 0.05)
%
% Other m-files required: pindrift, pinana, stabilize
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
% October 2017; Last revision: May 2018

%------------- BEGIN CODE --------------

p = asymstruct_eq.p;
asymstruct_newVapp = asymstruct_eq;

Vapp_array = 0:deltaV:Vend;
J_array = NaN(1, length(Vapp_array));

% first point is zero voltage
for i = 1:length(Vapp_array)
    [J_array(i), asymstruct_newVapp] = IgiveCurrentForVoltage(asymstruct_newVapp, Vapp_array(i));
    
    % eliminate negative values (e.g. the first value at zero Vapp can be negative even if it should be zero)
    if J_array(i) < 0
        disp([mfilename ' - Negative current of ' num2str(J_array(i)) ' mA/cm2 found for voltage ' num2str(Vapp_array(i)) ' V, setting to zero current'])
        J_array(i) = 0;
    end
end

LnJ_array = log(J_array);
prefactor = gradient(LnJ_array, Vapp_array);
nid_array = p.q ./ (p.kB .* p.T .* prefactor);

figure('Name', ['Ideality factor - ' inputname(1)], 'NumberTitle', 'off')
    plot(Vapp_array, nid_array, 'd-')
    xlabel('Applied voltage [V]');
    ylabel('Ideality factor');
    grid on
end

% this is adapted from findOptimVoc
function [current, asymstruct_newVapp] = IgiveCurrentForVoltage(asymstruct, Vapp)

    p = asymstruct.p;

    % no figures needed
    p.figson = 0;
    
    p.JV = 2; % mode for arbitrary Vapp profiles
    p.tpoints = 10;
    Vstart = p.Vapp; % current applied voltage
    Vend = Vapp; % new Voc    
    % the third parameter is the time point when the change from the old
    % voltage to the new one happens
    p.Vapp_params = [Vstart, Vend, 10 * p.t0];
    % starts at Vstart, linear Vapp decrease for the first points, then constant to Vend
    p.Vapp_func = @(coeff, t) (coeff(1) + (coeff(2)-coeff(1)).*t/coeff(3)) .* (t < coeff(3)) + coeff(2) .* (t >= coeff(3));
    
	asymstruct_newVapp = pindrift(asymstruct, p);

    % eliminate JV configuration before stabilizing
    p.JV = 0;
    % set Vapp as single value
    p.Vapp = Vend;

    asymstruct_newVapp = stabilize(asymstruct_newVapp);

    [~, Jn, ~] = pinana(asymstruct_newVapp);

    current = Jn(end);
    
end

%------------- END OF CODE --------------

