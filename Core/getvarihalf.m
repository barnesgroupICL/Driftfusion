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
varihalf = zeros(1, length(var)-1);
for i = 1:length(var)-1
    varihalf(1,i) = var(1,i)+((var(1,i+1)-var(1,i))/2);
end
end