function dev = build_device(par, meshoption)
% BUILD_DEVICE calls BUILD_PROPERTY for each device property. BUILD_PROPERTY then defines the
% properties at each point on the grid defined by MESHOPTION
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
switch meshoption
    case 'iwhole'
        xmesh = par.xx;
    case 'ihalf'
        xmesh = getvarihalf(par.xx);
end

% Constant properties
dev.taun = build_property(par.taun, xmesh, par, 'surface_rec_n', 0);
dev.taup = build_property(par.taup, xmesh, par, 'surface_rec_p', 0);

% Linearly graded properties
dev.EA = build_property(par.EA, xmesh, par, 'lin_graded', 0);
dev.IP = build_property(par.IP, xmesh, par, 'lin_graded', 0);
dev.mue = build_property(par.mue, xmesh, par, 'constant', 0);
dev.muh = build_property(par.muh, xmesh, par, 'constant', 0);
dev.mucat = build_property(par.mucat, xmesh, par, 'zeroed', 0);
dev.muani = build_property(par.muani, xmesh, par, 'zeroed', 0);
dev.epp = build_property(par.epp, xmesh, par, 'constant', 0);
dev.E0 = build_property(par.E0, xmesh, par, 'lin_graded', 0);
dev.Et = build_property(par.Et, xmesh, par, 'constant', 0);

% Logarithmically graded properties
dev.NA = build_property(par.NA, xmesh, par, 'log_graded', 0);
dev.ND = build_property(par.ND, xmesh, par, 'log_graded', 0);
dev.Nc = build_property(par.Nc, xmesh, par, 'log_graded', 0);
dev.Nv = build_property(par.Nv, xmesh, par, 'log_graded', 0);
dev.Nani = build_property(par.Nani, xmesh, par, 'log_graded', 0);
dev.Ncat = build_property(par.Ncat, xmesh, par, 'log_graded', 0);
dev.ni = build_property(par.ni, xmesh, par, 'log_graded', 0);
dev.n0 = build_property(par.n0, xmesh, par, 'log_graded', 0);
dev.p0 = build_property(par.p0, xmesh, par, 'log_graded', 0);
dev.DOSani = build_property(par.amax, xmesh, par, 'log_graded', 0);
dev.DOScat = build_property(par.cmax, xmesh, par, 'log_graded', 0);

% Properties that are zeroed in the interfaces
dev.B = build_property(par.B, xmesh, par, 'zeroed', 0);
dev.g0 = build_property(par.g0, xmesh, par, 'zeroed', 0);

% Gradient properties
dev.gradEA = build_property(par.EA, xmesh, par, 'lin_graded', 1);
dev.gradIP = build_property(par.IP, xmesh, par, 'lin_graded', 1);
dev.gradNc = build_property(par.Nc, xmesh, par, 'log_graded', 1);
dev.gradNv = build_property(par.Nv, xmesh, par, 'log_graded', 1);

dev.nt = distro_fun.nfun(dev.Nc, dev.EA, dev.Et, par.T, par.prob_distro_function);
dev.pt = distro_fun.pfun(dev.Nv, dev.IP, dev.Et, par.T, par.prob_distro_function);

end