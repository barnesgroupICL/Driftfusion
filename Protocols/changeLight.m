function struct_Int = changeLight(struct, newInt, tmax)
%CHANGELIGHT - Stabilize solutions at a new light intensity
%
% Syntax:  struct_Int = changeLight(struct, newInt, tmax)
%
% Inputs:
%   STRUCT - a solution struct as created by PINDRIFT.
%   NEWINT - the requested light intensity, zero is not supported as it is
%     much more robust to obtain a dark solution directly from equilibrate
%   TMAX - the initial stabilization time, can be zero for an automatic
%     guess
%
% Outputs:
%   STRUCT_INT - a solution struct at NEWINT light intensity
%
% Example:
%   ssol_i_1S_SR_1mS = changeLight(ssol_i_1S_SR, 1e-3, 5)
%     take the solution ssol_i_light and stabilize to a new light intensity
%     of 0.001, use as stabilization time 5 seconds
%   ssol_i_1S_SR_1mS = changeLight(ssol_i_1S_SR, 1e-3, 0)
%     as above, but estimate a good time for stabilization
%

% Other m-files required: df, verifyStabilization
% Subfunctions: none
% MAT-files required: none
%
% See also genIntStructs, df.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: January 2018

%------------- BEGIN CODE --------------

par = struct.par;
par.pulseon = 0;
par.tmesh_type = 2;
par.t0 = 1e-10;
par.tpoints = 30;
struct_Int = struct;

% set an initial time for stabilization tmax
if tmax
    tmax_temp = tmax;
else % if tmax was zero, estimate a good one
    if par.muion
        tmax_temp = min(1, 2^(-log10(par.muion)) / 10 + 2^(-log10(par.mue(1))));
    else
        tmax_temp = min(1e-3, 2^(-log10(par.mue(1))));
    end
end
par.tmax = tmax_temp;

% find the initial illumination intensity, if it's zero (dark) just use
% 1e-4, so that the first stabilization will be from dark to 1e-3
if par.Int
    oldInt = par.Int;
    steps = 1 + ceil(abs(log10(newInt / oldInt)));
else
    oldInt = 1e-3;
    steps = 1 + ceil(max(1, log10(newInt / oldInt)));
end

% if the step is just one, then newInt is the output of logspace function
Int_array = logspace(log10(oldInt), log10(newInt), steps);

%% change light in steps
% not needed to reach a good stabilized solution in
% each step, so stabilization is not verified here

% skip first value in the array as is the initial Int (and does not get through pindrift again) or, in case the input
% was in dark, the first value is 1e-3 and gets skipped
for i = 2:length(Int_array)
    disp([mfilename ' - Go from light intensity ' num2str(par.Int) ' to ' num2str(Int_array(i)) ' over ' num2str(par.tmax) ' s'])
    par.Int = Int_array(i); % set new light intensity
    struct_Int = df(struct_Int, par);
end

%% stabilize
struct_Int = stabilize(struct_Int); % go to steady state

%------------- END OF CODE --------------
