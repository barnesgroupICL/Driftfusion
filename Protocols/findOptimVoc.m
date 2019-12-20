function [sol_Voc, Voc] = findOptimVoc(sol)
%findOptimVOC - Stabilize an asymmetric solution to the real open circuit voltage
% which is the applied voltage where the current is minimized,
% starting from the voltage present in the given solution. The suggested
% workflow is to stabilize a symmetric solution (so that the charges
% profile are close to OC conditions) and break in a half with
% asymmetricize. Requires MATLAB's Optimization Toolbox.
%
% Syntax:  [sol_Voc, Voc] = findOptimVoc(sol)
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
% Other m-files required: df, pinAna, stabilize
% Subfunctions: IgiveCurrentForVoltage
% MAT-files required: none
%
% See also df, asymmetricize, findVoc, pinAna.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% April 2018; Last revision: April 2018

%------------- BEGIN CODE --------------

%assert(logical(sol.par.calcJ), [mfilename ' - calcJ needs to be set and different from zero']);

sol.par.tpoints = 10;
sol.par.JV = 0;

% find residual current
[~, originalCurrent] = dfana.calcJ(sol);

disp([mfilename ' - Original voltage: ' num2str(sol.par.Vapp, 8) ' V; original current: ' num2str(originalCurrent(end)) ' mA/cm2'])

if ~originalCurrent(end) % if the original current is exactly zero
    warning([mfilename ' - The initial residual current is zero, no better VOC can be found'])
    sol_Voc = sol;
    Voc = sol.par.Vapp;
    return;
end

% set an initial time for stabilization tmax
if max(sol.par.mucat)
    sol.par.tmax = min(5e0, 2^(-log10(max(sol.par.mucat))) / 10 + 2^(-log10(sol.par.mue(1))));
else
    sol.par.tmax = min(1e-2, 2^(-log10(sol.par.mue(1))));
end

sol.par.t0 = sol.par.tmax / 1e6;

%% estimate search range, assuming that a higher voltage causes a more positive current

% preallocate with the value from input solution
previousCurrent = originalCurrent(end);
% set a fallback value in case the for cycle doesn't find a good one
dVlimit = 1.3;
% look for the residual current at a new voltage,
% if the sign of the residual current is different use the voltage
% variation as a search range limit
for dV = [0.01, repelem(0.2, 4), repelem(0.05, 9)]
    % this assumes that a more positive voltage results in more positive
    % current
    [~, newCurrent, sol_newVapp] = IgiveCurrentForVoltage(sol, sol.par.Vapp - sign(previousCurrent) * dV);
    if sign(newCurrent) ~= sign(originalCurrent(end))
        dVlimit = dV;
        % if the current is smaller, the new solution is used as starting
        % point, otherwise the solution of the previous cycle is used
        if abs(newCurrent) < abs(previousCurrent)
            sol = sol_newVapp;
        end
        break
    else
        % the new solution is used as next point
        sol = sol_newVapp;
        % store the new solution current value for next cycle
        previousCurrent = newCurrent;
    end
end

disp([mfilename ' - Search range: ' num2str(dVlimit) ' V'])

%% do the search

% the tmax obtained in the previous step could be too big and make some next step to fail
sol.par.tmax = sol.par.tmax / 10;
sol.par.t0 = sol.par.t0 / 10;

% here sol is a fixed input, while Vapp is the input that will be
% optimized
fun = @(Vapp) IgiveCurrentForVoltage(sol, Vapp);

% StepTolerance is usually what stops the minimization here
% MaxFunctionEvaluations is 100 by default, reducing to 50 but even smaller
% values could be enough
options = optimoptions('fmincon', 'StepTolerance', 1e-7, 'FunctionTolerance', 1e-11, 'OptimalityTolerance', 1e-11, 'Algorithm', 'active-set', 'Display', 'notify', 'MaxFunctionEvaluations', 50);

% the starting point is the currently present voltage
% the constraints does not work when using the default interior-point
% algorithm, no idea why
Voc = fmincon(fun, sol.par.Vapp, [1; -1], [sol.par.Vapp + dVlimit; -(sol.par.Vapp - dVlimit)], [], [], [], [], [], options);

%% stabilize at the real Voc

sol.par.tpoints = 30;
sol.par.JV = 2; % mode for arbitrary Vapp profiles
Vstart = sol.par.Vapp; % current applied voltage
Vend = Voc; % new Voc
sol.par.tmax = sol.par.tmax / 10;
% the third parameter is the time point when the change from the old
% voltage to the new one happens
sol.par.Vapp_params = [Vstart, Vend, 10 * sol.par.t0];
% starts at Vstart, linear Vapp decrease for the first points, then constant to Vend
sol.par.Vapp_func = @(coeff, t) (coeff(1) + (coeff(2)-coeff(1)).*t/coeff(3)) .* (t < coeff(3)) + coeff(2) .* (t >= coeff(3));

sol_Voc = df(sol, sol.par);

%% stabilize

% eliminate JV configuration
sol_Voc.par.JV = 0;
% should be set anyway, but for avoiding risks, set Vapp
sol_Voc.par.Vapp = Vend;

sol_Voc = stabilize(sol_Voc); % go to steady state

% save Vend also in Vapp field of the struct
sol_Voc.Vapp = Vend;
% restore figson
sol_Voc.par.figson = 1;

disp([mfilename ' - VOC found at ' num2str(Vend, 8) ' V'])

end

%% take out the last current point
function [minimizeMe, current, sol_newVapp] = IgiveCurrentForVoltage(sol, Vapp)

par = sol.par;

par.JV = 2; % mode for arbitrary Vapp profiles
Vstart = par.Vapp; % current applied voltage
Vend = Vapp; % new Voc    
% the third parameter is the time point when the change from the old
% voltage to the new one happens
par.Vapp_params = [Vstart, Vend, 10 * par.t0];
% starts at Vstart, linear Vapp decrease for the first points, then constant to Vend
par.Vapp_func = @(coeff, t) (coeff(1) + (coeff(2)-coeff(1)).*t/coeff(3)) .* (t < coeff(3)) + coeff(2) .* (t >= coeff(3));

% the voltage step simulation can fail, shortening tmax (as in the catch block) could solve
try
    sol_newVapp = df(sol, par);
catch
    disp([mfilename ' - the voltage change failed, trying again with shorter tmax'])
    par.tmax = par.tmax / 10;
    par.t0 = par.t0 / 10;
    sol_newVapp = df(sol, par);
end

% eliminate JV configuration before stabilizing
par.JV = 0;
% set Vapp as single value
par.Vapp = Vend;

warning('off', 'df:verifyStabilization');
while ~verifyStabilization(sol_newVapp.sol, sol_newVapp.t, 1e-2) % check stability
    disp([mfilename ' - Stabilizing over ' num2str(par.tmax) ' s']);
    sol_newVapp = df(sol_newVapp, par);
    par.tmax = par.tmax * 5;
    par.t0 = par.t0 * 5;
end
warning('on', 'df:verifyStabilization');

[~, Jn] = dfana(sol_newVapp);

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
