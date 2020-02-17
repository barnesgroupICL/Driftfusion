function ppos = getpointpos(xpos, xmesh)
% Obtain point position PPOS from x-position XPOS for mesh defined by XMESH
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
if xpos <= xmesh(1)
    ppos = 1;
else
    ppos = find(xmesh <= xpos);
    ppos = ppos(end);
end

end
