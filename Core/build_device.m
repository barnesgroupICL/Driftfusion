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
dev.mucat = build_property(par.mucat, xmesh, par, 'constant', 0);
dev.muani = build_property(par.muani, xmesh, par, 'constant', 0);
dev.epp = build_property(par.epp, xmesh, par, 'constant', 0);
dev.nt = build_property(par.nt, xmesh, par, 'constant', 0);
dev.pt = build_property(par.pt, xmesh, par, 'constant', 0);
dev.sn = build_property(par.sn, xmesh, par, 'constant', 1);
dev.sp = build_property(par.sp, xmesh, par, 'constant', 1);
dev.ni = build_property(par.ni, xmesh, par, 'constant', 0);

% Linearly graded properties
dev.EA = build_property(par.EA, xmesh, par, 'lin_graded', 0);
dev.IP = build_property(par.IP, xmesh, par, 'lin_graded', 0);
dev.E0 = build_property(par.E0, xmesh, par, 'lin_graded', 0);

% Logarithmically graded properties
dev.Nc = build_property(par.Nc, xmesh, par, 'log_graded', 0);
dev.Nv = build_property(par.Nv, xmesh, par, 'log_graded', 0);
dev.n0 = build_property(par.n0, xmesh, par, 'log_graded', 0);
dev.p0 = build_property(par.p0, xmesh, par, 'log_graded', 0);
dev.Nani = build_property(par.Nani, xmesh, par, 'log_graded', 0);
dev.Ncat = build_property(par.Ncat, xmesh, par, 'log_graded', 0);
dev.DOSani = build_property(par.amax, xmesh, par, 'log_graded', 0);
dev.DOScat = build_property(par.cmax, xmesh, par, 'log_graded', 0);
dev.taun = build_property(par.taun, xmesh, par, 'log_graded', 0);
dev.taup = build_property(par.taup, xmesh, par, 'log_graded', 0);

% Properties that are zeroed in the interfaces
dev.g0 = build_property(par.g0, xmesh, par, 'zeroed', 0);
dev.B = build_property(par.B, xmesh, par, 'zeroed', 0);
dev.NA = build_property(par.NA, xmesh, par, 'zeroed', 0);
dev.ND = build_property(par.ND, xmesh, par, 'zeroed', 0);

% Gradient properties
dev.gradEA = build_property(par.EA, xmesh, par, 'lin_graded', 1);
dev.gradIP = build_property(par.IP, xmesh, par, 'lin_graded', 1);
dev.gradNc = build_property(par.Nc, xmesh, par, 'log_graded', 1);
dev.gradNv = build_property(par.Nv, xmesh, par, 'log_graded', 1);

% Surface recombination velocity equivalence schemes
%dev.mue = build_property(par.mue, xmesh, par, 'mue_interface', 0);
%dev.muh = build_property(par.muh, xmesh, par, 'muh_interface', 0);
dev.mue = build_property(par.mue, xmesh, par, 'constant', 0);
dev.muh = build_property(par.muh, xmesh, par, 'constant', 0);

dev.taun_vsr = build_property(par.taun, xmesh, par, 'taun_vsr', 0);
dev.taup_vsr = build_property(par.taup, xmesh, par, 'taup_vsr', 0);
dev.alpha_prime = build_property(par.xx, xmesh, par, 'alpha_prime', 1);
dev.beta_prime = build_property(par.xx, xmesh, par, 'beta_prime', 1);

% Switches
dev.int_switch = build_property(par.int_switch, xmesh, par, 'int_switch', 1);
dev.bulk_switch = abs(dev.int_switch-1);

% Tranlsated co-ordinates
dev.xprime = build_property(par.xx, xmesh, par, 'xprime', 1);
dev.xprime_n = build_property(par.xx, xmesh, par, 'xprime_n', 1);
dev.xprime_p = build_property(par.xx, xmesh, par, 'xprime_p', 1);
end
