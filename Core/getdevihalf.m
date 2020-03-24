function devihalf = getdevihalf(par)
% Builds the device structure on the subinterval XMESH
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
dev = par.dev;
fnames = fieldnames(dev);

for i = 1:length(fnames)
    prop = char(fnames(i));
    eval(['devihalf.',prop,' = dev.',prop,'(1:end-1) + (dev.',prop,'(2:end) - dev.',prop,'(1:end-1))/2;']);
end



end