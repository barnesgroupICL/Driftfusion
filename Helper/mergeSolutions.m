function sol = mergeSolutions(varargin)
%MERGESOLUTIONS - Concatenates solutions along the time dimension
% 
% Syntax:  sol = mergeSolutions(varargin)
%
% Inputs:
%   VARARGIN - various solutions structures as created by DF
%
% Outputs:
%   SOL - the solution obtained joining all the provided solutions over the time dimension
%
% Example:
%   sol = mergeSolutions(soleq.ion, sol_ramp, sol_dwell)
%     joins the three solutions so that the output will show the temporal
%     evolution of all of these
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also df.
%
%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%------------- BEGIN CODE --------------

xPointsNumber = zeros(1, length(varargin));
varsNumber = zeros(1, length(varargin));
for i = 1:length(varargin)
    xPointsNumber(i) = size(varargin{i}.u, 2);
    varsNumber(i) = size(varargin{i}.u, 3);
end
assert(~range(xPointsNumber), [mfilename ' - The provided solutions does not have the same number of spatial points.']);
assert(~range(varsNumber), [mfilename ' - The provided solutions does not have the same number of variables.']);

% use the first time point of the first input for preallocating
sol.u = varargin{i}.u(1,:,:);
sol.t = varargin{i}.t(1);
tpoints = 1;

sol.x = varargin{end}.x;
sol.par = varargin{end}.par;
sol.par.t0 = varargin{1}.par.t0;
start_time = 0;

for i = 1:length(varargin)
    % discard all the first timepoints for each input: the one from the first
    % input is already included during the preallocation, while the others
    % can be discarded as they (should be) just duplicate of the last
    % timepoint of the previous input
    sol.u = cat(1, sol.u, varargin{i}.u(2:end,:,:));
    tmesh = start_time + varargin{i}.t(2:end);
    sol.t = [sol.t, tmesh];
    start_time = tmesh(end);
    % subtract one as we're excluding the first timepoint
    tpoints = tpoints + varargin{i}.par.tpoints - 1;
end
sol.par.tmax = start_time;
sol.par.tpoints = tpoints;

end

%------------- END OF CODE --------------