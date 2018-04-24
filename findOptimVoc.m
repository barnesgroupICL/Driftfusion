function [asymstruct_voc, Voc] = findOptimVoc(asymstruct)
%findOptimVOC - Stabilize an asymmetrical solution quite at Voc to the real open circuit voltage
% which is the applied voltage where the current is minimized,
% starting from the voltage present in the given solution and scanning in a
% 20 mV range. Needs Optimization Toolbox.
%
% Syntax:  [asymstruct_voc, Voc] = findOptimVoc(asymstruct)
%
% Inputs:
%   ASYMSTRUCT - an asymmetric struct as created by PINDRIFT using the OC
%     parameter set to 0, ideally breaking a symmetric solution using
%     asymmetricize function so that the starting Vapp is aready close to
%     the real VOC.
%
% Outputs:
%   ASYMSTRUCT_VOC - an asymmetric struct as created by PINDRIFT using the OC
%     parameter set to 0 and applied voltage identical to the real open
%     circuit voltage
%
% Example:
%   sol_i_light_SR_OC = findOptimVoc(asymmetricize(ssol_i_light_SR, 1))
%     get closer to the real VOC
%
% Other m-files required: pindrift
% Subfunctions: IgiveCurrentForVoltage
% MAT-files required: none
%
% See also pindrift.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% April 2018; Last revision: April 2018

%------------- BEGIN CODE --------------

asymstruct.p.figson = 0;
asymstruct.p.tpoints = 10;
asymstruct.p.Ana = 1;
asymstruct.p.calcJ = 1;

% set an initial time for stabilization tmax
if asymstruct.p.mui
    asymstruct.p.tmax = min(5e1, 2^(-log10(asymstruct.p.mui)) / 10 + 2^(-log10(asymstruct.p.mue_i)));
else
    asymstruct.p.tmax = min(5e1, 2^(-log10(asymstruct.p.mue_i)));
end

asymstruct.p.t0 = asymstruct.p.tmax / 1e6;

% here asymstruct is a fixed input, while Vapp is the input that will be
% optimized
fun = @(Vapp) IgiveCurrentForVoltage(asymstruct, Vapp);

% StepTolerance is usually what stops the minimization here
options = optimoptions('fmincon', 'StepTolerance', 1e-7, 'FunctionTolerance', 1e-13, 'OptimalityTolerance', 1e-13, 'Algorithm', 'active-set');

% the starting point is the currently present voltage
% the constraints does not work when using the default interior-point
% algorithm, no idea why
Voc = fmincon(fun, asymstruct.p.Vapp, [1; -1], [asymstruct.p.Vapp + 0.01; -(asymstruct.p.Vapp - 0.01)], [], [], [], [], [], options);

%% stabilize at the real Voc

asymstruct.p.tpoints = 30;
asymstruct.p.JV = 2; % mode for arbitrary Vapp profiles
Vstart = asymstruct.p.Vapp; % current applied voltage
Vend = Voc; % new Voc    
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


warning('off', 'pindrift:verifyStabilization');
while ~verifyStabilization(asymstruct_voc.sol, asymstruct_voc.t, 1e-3) % check stability
    disp([mfilename ' - Stabilizing over ' num2str(asymstruct_voc.p.tmax) ' s']);
    asymstruct_voc = pindrift(asymstruct_voc, asymstruct_voc.p);
    asymstruct_voc.p.tmax = asymstruct_voc.p.tmax * 10;
end
warning('on', 'pindrift:verifyStabilization');

% save Vend also in Vapp field of the struct
asymstruct_voc.Vapp = Vend;
% restore figson
asymstruct_voc.p.figson = 1;

end

%% take out the last current point
function minimizeMe = IgiveCurrentForVoltage(asymstruct, Vapp)

p = asymstruct.p;

p.JV = 2; % mode for arbitrary Vapp profiles
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

warning('off', 'pindrift:verifyStabilization');
while ~verifyStabilization(asymstruct_newVapp.sol, asymstruct_newVapp.t, 1e-2) % check stability
    disp([mfilename ' - Stabilizing over ' num2str(p.tmax) ' s']);
    asymstruct_newVapp = pindrift(asymstruct_newVapp, p);
    p.tmax = p.tmax * 10;
end
warning('on', 'pindrift:verifyStabilization');

[~, ~, Jn, ~, ~, ~] = pinAna(asymstruct_newVapp);

% I want the absolute value, but it's nicer to take the squared rather than
% the abs, for not introducing non derivable points
% the stupid multiplicative factor is needed for getting far from the eps
% (double precision) which is the minimum FunctionTolerance which can be
% set
minimizeMe = 1e30 * Jn(end)^2;

disp([mfilename ' - Using with Vapp ' num2str(Vapp, 8) ' V, a current of ' num2str(Jn(end)) ' mA/cm2 was found'])

end

%------------- END OF CODE --------------