function xsolver = getxihalf(sol)
% Builds subinterval xmesh from xmesh- can use GETVARIHALF instead
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
x = sol.par.xx;
for i = 1:length(x)-1
    xsolver(i) = x(i)+((x(i+1)-x(i))/2);
end
end