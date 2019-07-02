function symstruct_pulse = TPV_single_exec(symstruct_Int, tmax, pulse_int)
%TPV_SINGLE_EXEC - Perform a single Transient PhotoVoltage (TPV) simulation
%
% Syntax:  symstruct_pulse = TPV_single_exec(symstruct_Int, t_max, pulse_int)
%
% Inputs:
%   SYMSTRUCT_INT - a single symmetric struct as created by PINDRIFT.
%   TMAX - a first guess for maximum time which should be enough for
%     simulating the full decay, if zero is provided a value of 1e-5 s is
%     used and no further exploration of tmax is performed (this is useful
%     when we are just interested in the amount of voltage variation,
%     rather than in the full decay)
%   PULSE_INT - intensity of the perturbing light pulse
%
% Outputs:
%   SYMSTRUCT_PULSE - a struct with a solution being perturbed by a light
%     pulse
%
% Example:
%   TPV_single_exec(ssol_i_1S_SR, 1e-3, 5)
%     simulate a TPV on the stabilized solution provided, a time window
%     of at least a millisecond, and a light pulse intensity of 5
%
% Other m-files required: pindrift, verifyStabilization
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift, TPVvariab_full_exec, TPVconst_full_exec, verifyStabilization, TPV_single_analysis.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

p = symstruct_Int.p;

if tmax % if tmax is provided different from zero
    p.tmax = tmax;
else % if tmax is set as zero means that reacing a good tmax is not required
    p.tmax = 1e-5;
end

p.tmesh_type = 3; % specific for TPV
p.tpoints = 200;

% set pulse characteristics
p.pulseon = 1;
p.pulselen = 1e-8; % 10 ns pulse length
p.pulseint = pulse_int;
p.pulsestart = p.pulselen;

p.t0 = min(p.tmax/1e5, p.pulselen/100); % needs to be smaller than pulselen

% pre-allocate variable to big number
lastIterationFinalVoltage = Inf;

%% do TPV
warning('off', 'pindrift:verifyStabilization'); % in this part the stability warning can be ignored
tmax_ok = false;
while ~tmax_ok
    disp([mfilename ' - TPV tmax: ' num2str(p.tmax) ' s']);
    % each cycle starts with the input solution which should be provided stabilized, with no pulse, so that light intensity is already ok for sure
    symstruct_pulse = pindrift(symstruct_Int, p);
    % verify if the decay reached a stable condition, and that the resulting
    % ending voltage is small

    % at low light intensities the solution could behave weirdly (can oscillate) and
    % verifyStabilization could fail.
    if symstruct_pulse.p.Int < 1e-6
        if symstruct_pulse.Voc(end) < symstruct_pulse.Voc(1)
            condition1 = true;
            disp('Weird things happening, final voltage lower than initial voltage, assuming that the cell is stabilized without checking')
        else
            condition1 = verifyStabilization(symstruct_pulse.sol, symstruct_pulse.t, 0.5);
        end
    else
        condition1 = verifyStabilization(symstruct_pulse.sol, symstruct_pulse.t, 0.05);
    end
    
    % this condition works well at high light intensity, but fails where
    % noise gets important (at low light intensities)
    if symstruct_pulse.p.Int > 1e-5
        threshold = 0.01;
    elseif symstruct_pulse.p.Int > 1e-7
        threshold = 0.03;
    else
        threshold = 0.1;
    end
    condition2 = (symstruct_pulse.Voc(end) - symstruct_pulse.Voc(1)) < threshold * (max(symstruct_pulse.Voc) - symstruct_pulse.Voc(1));
    if condition1 && condition2
        tmax_ok = true;
    elseif ~tmax
        % if tmax is set as zero means that reacing a good tmax is not required, for example when looking for a good pulse intensity
        tmax_ok = true;
    elseif symstruct_pulse.Voc(end) > lastIterationFinalVoltage
        % if the final voltage is higher than the previous iteration,
        % something wrong is happening (the longer the stabilization time
        % the lower the final voltage we would expect) so things are not
        % going to get better with even longer times, better stop here
        disp('Weird things happening, last iteration final voltage higher than previous iteration one, assuming that the cell is stabilized without checking')
        tmax_ok = true;
    else
        p.tmax = p.tmax * 5; % the output solution will include the last used value, not this one, so the actually used one
    end
    lastIterationFinalVoltage = symstruct_pulse.Voc(end);
end
warning('on', 'pindrift:verifyStabilization');

% verify that the pulse atually did something
illuminated_tpoints = symstruct_pulse.t(symstruct_pulse.t < p.pulselen);
if isempty(illuminated_tpoints) % this is really unlikely, as t0 is set to a value smaller than pulselen
    warning('The laser pulse was too short or the time mesh too loose');
end

%------------- END OF CODE --------------
