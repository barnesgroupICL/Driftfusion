function varihalf = getvarihalf(var)
% Builds individual variable arrays on the subinterval xmesh
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Start code

assert(isnumeric(var) && size(var,2) > 1, [mfilename ' - The input ' inputname(1)...
    ' has to be provided as a row vector or as a matrix'])

varihalf = (var(:,1:(end-1)) + var(:,2:end)) / 2;
end