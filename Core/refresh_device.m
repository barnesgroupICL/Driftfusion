function par = refresh_device(par)
% Rebuilds important device properties
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
par.xx = meshgen_x(par);
par.x_sub = getvar_sub(par.xx);
par.dev = build_device(par, 'whole');
par.dev_sub = build_device(par, 'sub');
% Get generation profiles
par.gx1 = generation(par, par.light_source1, par.laser_lambda1);
par.gx2 = generation(par, par.light_source2, par.laser_lambda2);

end