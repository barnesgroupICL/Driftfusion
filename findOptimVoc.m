function [asymstruct_voc, Voc] = findOptimVoc(asymstruct)
%findOptimVOC - Stabilize an asymmetric solution to the real open circuit voltage
% which is the applied voltage where the current is minimized,
% starting from the voltage present in the given solution. The suggested
% workflow is to stabilize a symmetric solution (so that the charges
% profile are close to OC conditions) and break in a half with
% asymmetricize. Requires MATLAB's Optimization Toolbox.
%
% Syntax:  [asymstruct_voc, Voc] = findOptimVoc(asymstruct)
%
% Inputs:
%   ASYMSTRUCT - an asymmetric struct as created by PINDRIFT, ideally
%     breaking a symmetric solution using asymmetricize function so that
%     the starting Vapp is aready close to the real VOC.
%
% Outputs:
%   ASYMSTRUCT_VOC - an asymmetric struct as created by PINDRIFT using the OC
%     parameter set to 0 and applied voltage identical to the real open
%     circuit voltage
%   VOC - the value of the obtained open circuit voltage
%
% Example:
%   [sol_i_1S_SR_OC, VOC_1S] = findOptimVoc(asymmetricize(ssol_i_1S_SR))
%     get closer to the real VOC and save the VOC value in VOC_1S variable
%
% Other m-files required: pindrift, pinAna, stabilize
% Subfunctions: IgiveCurrentForVoltage
% MAT-files required: none
%
% See also pindrift, asymmetricize, findVoc, pinAna.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% April 2018; Last revision: April 2018

%------------- BEGIN CODE --------------

assert(logical(asymstruct.p.calcJ), [mfilename ' - calcJ needs to be set and different from zero']);

asymstruct.p.figson = 0;
asymstruct.p.tpoints = 10;
asymstruct.p.Ana = 1;
asymstruct.p.JV = 0;

% find residual current
[~, ~, originalCurrent] = pinana(asymstruct);

disp([mfilename ' - Original voltage: ' num2str(asymstruct.p.Vapp, 8) ' V; original current: ' num2str(originalCurrent(end)) ' mA/cm2'])

if ~originalCurrent(end) % if the original current is exactly zero
    warning([mfilename ' - The initial residual current is zero, no better VOC can be found'])
    asymstruct_voc = asymstruct;
    Voc = asymstruct.p.Vapp;
    return;
end

% set an initial time for stabilization tmax
if max(asymstruct.p.muion)
    asymstruct.p.tmax = min(5e0, 2^(-log10(max(asymstruct.p.muion))) / 10 + 2^(-log10(asymstruct.p.mue(1))));
else
    asymstruct.p.tmax = min(1e-2, 2^(-log10(asymstruct.p.mue(1))));
end

asymstruct.p.t0 = asymstruct.p.tmax / 1e6;

%% estimate search range, assuming that a higher voltage causes a more positive current

% preallocate with the value from input solution
previousCurrent = originalCurrent(end);
% set a fallback value in case the for cycle doesn't find a good one
dVlimit = 1.3;
% look for the residual current at a new voltage,
% if the sign of the residual current is different use the voltage
% variation as a search range limit
for dV = [0.01, 0.7, repelem(0.05, 11)]
    % this assumes that a more positive voltage results in more positive
    % current
    [~, newCurrent, asymstruct_newVapp] = IgiveCurrentForVoltage(asymstruct, asymstruct.p.Vapp - sign(previousCurrent) * dV);
    if sign(newCurrent) ~= sign(originalCurrent(end))
        dVlimit = dV;
        % if the current is smaller, the new solution is used as starting
        % point, otherwise the solution of the previous cycle is used
        if abs(newCurrent) < abs(previousCurrent)
            asymstruct = asymstruct_newVapp;
        end
        break
    else
        % the new solution is used as next point
        asymstruct = asymstruct_newVapp;
        % store the new solution current value for next cycle
        previousCurrent = newCurrent;
    end
end

disp([mfilename ' - Search range: ' num2str(dVlimit) ' V'])

%% do the search

% the tmax obtained in the previous step could be too big and make some next step to fail
asymstruct.p.tmax = asymstruct.p.tmax / 10;
asymstruct.p.t0 = asymstruct.p.t0 / 10;

% here asymstruct is a fixed input, while Vapp is the input that will be
% optimized
fun = @(Vapp) IgiveCurrentForVoltage(asymstruct, Vapp);

% StepTolerance is usually what stops the minimization here
options = optimoptions('fmincon', 'StepTolerance', 1e-7, 'FunctionTolerance', 1e-13, 'OptimalityTolerance', 1e-13, 'Algorithm', 'active-set', 'Display', 'notify');

% the starting point is the currently present voltage
% the constraints does not work when using the default interior-point
% algorithm, no idea why
Voc = fmincon(fun, asymstruct.p.Vapp, [1; -1], [asymstruct.p.Vapp + dVlimit; -(asymstruct.p.Vapp - dVlimit)], [], [], [], [], [], options);

%% stabilize at the real Voc

asymstruct.p.tpoints = 30;
asymstruct.p.JV = 2; % mode for arbitrary Vapp profiles
Vstart = asymstruct.p.Vapp; % current applied voltage
Vend = Voc; % new Voc
asymstruct.p.tmax = asymstruct.p.tmax / 10;
% the third parameter is the time point when the change from the old
% voltage to the new one happens
asymstruct.p.Vapp_params = [Vstart, Vend, 10 * asymstruct.p.t0];
% starts at Vstart, linear Vapp decrease for the first points, then constant to Vend
asymstruct.p.Vapp_func = @(coeff, t) (coeff(1) + (coeff(2)-coeff(1)).*t/coeff(3)) .* (t < coeff(3)) + coeff(2) .* (t >= coeff(3));

asymstruct_voc = pindrift(asymstruct, asymstruct.p);

%% stabilize

% eliminate JV configuration
asymstruct_voc.p.JV = 0;
% should be set anyway, but for avoiding risks, set Vapp
asymstruct_voc.p.Vapp = Vend;

asymstruct_voc = stabilize(asymstruct_voc); % go to steady state

% save Vend also in Vapp field of the struct
asymstruct_voc.Vapp = Vend;
% restore figson
asymstruct_voc.p.figson = 1;

disp([mfilename ' - VOC found at ' num2str(Vend, 8) ' V'])

end

%% take out the last current point
function [minimizeMe, current, asymstruct_newVapp] = IgiveCurrentForVoltage(asymstruct, Vapp)

p = asymstruct.p;

p.JV = 2; % mode for arbitrary Vapp profiles
Vstart = p.Vapp; % current applied voltage
Vend = Vapp; % new Voc    
% the third parameter is the time point when the change from the old
% voltage to the new one happens
p.Vapp_params = [Vstart, Vend, 10 * p.t0];
% starts at Vstart, linear Vapp decrease for the first points, then constant to Vend
p.Vapp_func = @(coeff, t) (coeff(1) + (coeff(2)-coeff(1)).*t/coeff(3)) .* (t < coeff(3)) + coeff(2) .* (t >= coeff(3));

% the voltage step simulation can fail, shortening tmax (as in the catch block) could solve
try
    asymstruct_newVapp = pindrift(asymstruct, p);
catch
    disp([mfilename ' - the voltage change failed, trying again with shorter tmax'])
    p.tmax = p.tmax / 10;
    p.t0 = p.t0 / 10;
    asymstruct_newVapp = pindrift(asymstruct, p);
end

% eliminate JV configuration before stabilizing
p.JV = 0;
% set Vapp as single value
p.Vapp = Vend;

warning('off', 'pindrift:verifyStabilization');
while ~verifyStabilization(asymstruct_newVapp.sol, asymstruct_newVapp.t, 1e-2) % check stability
    disp([mfilename ' - Stabilizing over ' num2str(p.tmax) ' s']);
    asymstruct_newVapp = pindrift(asymstruct_newVapp, p);
    p.tmax = p.tmax * 5;
    p.t0 = p.t0 * 5;
end
warning('on', 'pindrift:verifyStabilization');

[~, ~, Jn] = pinana(asymstruct_newVapp);

current = Jn(end);
% I want the absolute value, but it's nicer to take the squared rather than
% the abs, for not introducing non derivable points
% the stupid multiplicative factor is needed for getting far from the eps
% (double precision) which is the minimum FunctionTolerance which can be
% set
minimizeMe = 1e30 * current^2;

disp([mfilename ' - Using with Vapp ' num2str(Vapp, 8) ' V, a current of ' num2str(Jn(end)) ' mA/cm2 was found'])

end

%------------- END OF CODE --------------
