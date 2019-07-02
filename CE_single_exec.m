function asymstruct_CE = CE_single_exec(asymstruct_Int, tmax)
%CE_SINGLE_EXEC - Do a single Charge Extraction (CE) experiment
%
% Syntax:  asymstruct_CE = CE_single_exec(asymstruct_Int, t_max)
%
% Inputs:
%   ASYMSTRUCT_INT - a single asymmetric struct as created by PINDRIFT.
%   TMAX - a first guess for maximum time which should be enough for
%     registering the full extraction
%
% Outputs:
%   ASYMSTRUCT_CE - a struct with a solution evolving after the start of
%     the CE experiment
%
% Example:
%   CE_single_exec(asymmetricize(ssol_i_1S_SR), 1e-3)
%     simulate CE on the stabilized solution provided, a time window
%     of at least a millisecond, and blocking contacts boundary conditions
%
% Other m-files required: pindrift, verifyStabilization
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift, CE_full_exec, verifyStabilization.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

p = asymstruct_Int.p;

p.tmesh_type = 2;
p.tmax = tmax;
p.t0 = 1e-10;
p.tpoints = 20; % rough, just for finding tmax
p.calcJ = 0; % we need the current for calculating the charge extracted, but we calculate it just when we got the right tmax
p.Int = 0; % the current extraction happens in dark

% ideal sudden short circuit cannot be solved
% boundary conditions have to be changed more slowly 
% the solver needs at least one point at the starting voltage, before switching to short circuit
p.JV = 2; % mode for arbitrary Vapp changes
p.Vstart = asymstruct_Int.p.Vapp;
p.Vend = 0; % final Voc, short circuit

% the third parameter is the time point when the change from the old
% voltage to the new one happens
p.Vapp_params = [p.Vstart, p.Vend, 10 * asymstruct_Int.p.t0];
% starts at Vstart, linear Vapp decrease for the first points, then constant to Vend
p.Vapp_func = @(coeff, t) (coeff(1) + (coeff(2)-coeff(1)).*t/coeff(3)) .* (t < coeff(3)) + coeff(2) .* (t >= coeff(3));

%% look for a long enough tmax
warning('off','pindrift:verifyStabilization'); % in this part the stability warning can be ignored

t_max_ok = false;
while ~t_max_ok
    disp([mfilename ' - Int pre CE: ' num2str(asymstruct_Int.p.Int) '; tmax: ' num2str(p.tmax) ' s']);
    % each cycle starts with the input solution which should be provided stabilized, with no pulse, so that light intensity is already ok for sure
    asymstruct_CE = pindrift(asymstruct_Int, p);
    % verify if the decay reached a stable condition
    if verifyStabilization(asymstruct_CE.sol, asymstruct_CE.t, 0.25)
        t_max_ok = true;
    else
        p.tmax = p.tmax * 10; % the solution will include the last used value, so the actually used one
    end
end

%% do the actual CE
p.tpoints = 200;

% we need the current to be included in the solution
p.Ana = 1;
if ~p.calcJ
    p.calcJ = 2;
end

asymstruct_CE = pindrift(asymstruct_Int, p);

warning('on','pindrift:verifyStabilization'); % re-enable the warnings

%------------- END OF CODE --------------
