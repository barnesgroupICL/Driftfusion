function [struct_voc, Voc] = findVocOptim(struct, approxVoc)
%FINDVOCOPTIM - Stabilize a solution to the real open circuit voltage using MATLAB's Optimization Toolbox
% The Voc found is the applied voltage where the current is minimized,
% starting from the voltage present in the given solution.
% Requires MATLAB's Optimization Toolbox.
%
% Syntax:  [struct_voc, Voc] = findVocOptim(struct)
%
% Inputs:
%   STRUCT - a struct as created by DF
%   APPROXVOC - optional numeric, an initial guess for the Voc voltage, if
%     it is missing, the search will take more time
%
% Outputs:
%   STRUCT_VOC - a struct as created by DF with applied voltage identical
%     to the real open circuit voltage
%   VOC - the value of the obtained open circuit voltage
%
% Example:
%   [soleq_ion_voc, Voc] = findVocOptim(soleq.ion)
%     get closer to the real VOC and save the open circuit voltage value in Voc variable
%
% Other m-files required: df, dfana, stabilize
% Subfunctions: IgiveCurrentForVoltage
% MAT-files required: none
%
% See also df, findVoc, findVocDirect, dfana.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% April 2018; Last revision: June 2020

%------------- BEGIN CODE --------------


% set an initial time for stabilization tmax
if any(struct.par.mu_c) && struct.par.mobseti
    struct.par.tmax = min(5e0, 2^(-log10(max(struct.par.mu_c))) / 10 + 2^(-log10(min(struct.par.mu_n))));
else
    struct.par.tmax = min(1e-2, 2^(-log10(min(struct.par.mu_n))));
end

struct.par.tmesh_type = 2;
struct.par.tpoints = 10;
struct.par.t0 = struct.par.tmax / 1e6;

% if an approximated Voc is provided, start from that
if nargin > 1 && getVend(struct) ~= approxVoc
    struct.par.V_fun_arg = [getVend(struct), approxVoc, 10 * struct.par.t0];
    % starts at Vstart, linear Vapp decrease for the first points, then constant to Vend
    struct.par.V_fun_type = 'sweepAndStill';
    struct = df(struct, struct.par);

    struct.par.V_fun_type = 'constant';
    struct.par.V_fun_arg = approxVoc;
    struct = stabilize(struct); % go to steady state
end

% find residual current
originalCurrent = dfana.calcJ(struct);

disp([mfilename ' - Original voltage: ' num2str(getVend(struct), 8) ' V; original current: ' num2str(originalCurrent.tot(end,end)) ' A/cm2'])

%% estimate search range, assuming that a higher voltage causes a more positive current

% preallocate with the value from input solution
previousCurrent = originalCurrent.tot(end,end);
% set a fallback value in case the for cycle doesn't find a good one
dVlimit = 1.4;
% look for the residual current at a new voltage,
% if the sign of the residual current is different use the voltage
% variation as a search range limit
for dV = [0.001, 0.01, 0.05, repelem(0.2, 4), repelem(0.05, 9)]
    % this assumes that a more positive voltage results in more positive
    % current
    [~, newCurrent, struct_newVapp] = IgiveCurrentForVoltage(struct, struct.par.Vapp - sign(previousCurrent) * dV);
    if sign(newCurrent) ~= sign(originalCurrent.tot(end,end))
        dVlimit = dV;
        % if the current is smaller, the new solution is used as starting
        % point, otherwise the solution of the previous cycle is used
        if abs(newCurrent) < abs(previousCurrent)
            struct = struct_newVapp;
        end
        break
    else
        % the new solution is used as next point
        struct = struct_newVapp;
        % store the new solution current value for next cycle
        previousCurrent = newCurrent;
    end
end

disp([mfilename ' - Search range: ' num2str(dVlimit) ' V'])

%% do the search

% the tmax obtained in the previous step could be too big and make some next step to fail
struct.par.tmax = struct.par.tmax / 10;
struct.par.t0 = struct.par.t0 / 10;

% here struct is a fixed input, while Vapp is the input that will be
% optimized
fun = @(Vapp) IgiveCurrentForVoltage(struct, Vapp);

% StepTolerance is usually what stops the minimization here
% MaxFunctionEvaluations is 100 by default, reducing to 30 but even smaller
% values could be enough
options = optimoptions(@fmincon, 'StepTolerance', 1e-7, 'FunctionTolerance', 1e-11, 'OptimalityTolerance', 1e-11, 'Algorithm', 'active-set', 'Display', 'notify', 'MaxFunctionEvaluations', 30);

% the starting point is the currently present voltage
% the constraints does not work when using the default interior-point
% algorithm, no idea why
Voc = fmincon(fun, getVend(struct), [1; -1], [struct.par.Vapp + dVlimit; -(struct.par.Vapp - dVlimit)], [], [], [], [], [], options);

%% stabilize at the real Voc

struct.par.tpoints = 20;

struct.par.V_fun_arg = [getVend(struct), Voc, 10 * struct.par.t0];
% starts at Vstart, linear Vapp decrease for the first points, then constant to Vend
struct.par.V_fun_type = 'sweepAndStill';

struct_voc = df(struct, struct.par);

%% stabilize

struct_voc.par.V_fun_type = 'constant';
struct_voc.par.V_fun_arg = Voc;

struct_voc = stabilize(struct_voc); % go to steady state

residualJ = dfana.calcJ(struct_voc);

disp([mfilename ' - VOC found at ' num2str(Voc, 10) ' V, residual current ' num2str(residualJ.tot(end,end)) ' A/cm2'])

end

%% take out the last current point
function [minimizeMe, current, struct_newVapp] = IgiveCurrentForVoltage(struct, Vapp)

par = struct.par;

Vstart = par.Vapp; % current applied voltage
Vend = Vapp; % new Voc
% the third parameter is the time point when the change from the old
% voltage to the new one happens
par.V_fun_arg = [Vstart, Vend, 10 * par.t0];
% starts at Vstart, linear Vapp decrease for the first points, then constant to Vend
par.V_fun_type = 'sweepAndStill';

% the voltage step simulation can fail, shortening tmax (as in the catch block) could solve
try
    struct_newVapp = df(struct, par);
catch
    disp([mfilename ' - the voltage change failed, trying again with shorter tmax'])
    par.tmax = par.tmax / 10;
    par.t0 = par.t0 / 10;
    struct_newVapp = df(struct, par);
end

struct_newVapp.par.V_fun_type = 'constant';
struct_newVapp.par.V_fun_arg = Vend;

struct_newVapp = stabilize(struct_newVapp);

J = dfana.calcJ(struct_newVapp);

current = J.tot(end,end);
% I want the absolute value, but it's nicer to take the squared rather than
% the abs, for not introducing non derivable points
% the stupid multiplicative factor is needed for getting far from the eps
% (double precision) which is the minimum FunctionTolerance which can be
% set
minimizeMe = 1e30 * current^2;

disp([mfilename ' - Using with Vapp ' num2str(getVend(struct_newVapp), 8) ' V, a current of ' num2str(current) ' A/cm2 was found'])

end

%------------- END OF CODE --------------
