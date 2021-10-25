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
    case 'whole'
        xmesh = par.xx;
    case 'sub'
        xmesh = par.x_sub;
end

% Constant properties
dev.mu_c = build_property(par.mu_c, xmesh, par, 'constant', 0);
dev.mu_a = build_property(par.mu_a, xmesh, par, 'constant', 0);

dev.sn = build_property(par.sn, xmesh, par, 'constant', 1);
dev.sp = build_property(par.sp, xmesh, par, 'constant', 1);
dev.mu_n = build_property(par.mu_n, xmesh, par, 'constant', 0);
dev.mu_p = build_property(par.mu_p, xmesh, par, 'constant', 0);
    
% Linearly graded properties
dev.Phi_EA = build_property(par.Phi_EA, xmesh, par, 'lin_graded', 0);
dev.Phi_IP = build_property(par.Phi_IP, xmesh, par, 'lin_graded', 0);
dev.EF0 = build_property(par.EF0, xmesh, par, 'lin_graded', 0);

% Exponentially graded properties
dev.Nc = build_property(par.Nc, xmesh, par, 'exp_graded', 0);
dev.Nv = build_property(par.Nv, xmesh, par, 'exp_graded', 0);
dev.n0 = build_property(par.n0, xmesh, par, 'exp_graded', 0);
dev.p0 = build_property(par.p0, xmesh, par, 'exp_graded', 0);
dev.Nani = build_property(par.Nani, xmesh, par, 'exp_graded', 0);
dev.Ncat = build_property(par.Ncat, xmesh, par, 'exp_graded', 0);
dev.a_max = build_property(par.a_max, xmesh, par, 'exp_graded', 0);
dev.c_max = build_property(par.c_max, xmesh, par, 'exp_graded', 0);
dev.taun = build_property(par.taun, xmesh, par, 'exp_graded', 0);
dev.taup = build_property(par.taup, xmesh, par, 'exp_graded', 0);

% Properties that are zeroed in the interfaces
dev.g0 = build_property(par.g0, xmesh, par, 'zeroed', 0);
dev.B = build_property(par.B, xmesh, par, 'zeroed', 0);
dev.NA = build_property(par.NA, xmesh, par, 'exp_graded', 0);
dev.ND = build_property(par.ND, xmesh, par, 'exp_graded', 0);

% Gradient properties
dev.gradEA = build_property(par.Phi_EA, xmesh, par, 'lin_graded', 1);
dev.gradIP = build_property(par.Phi_IP, xmesh, par, 'lin_graded', 1);
dev.gradNc = build_property(par.Nc, xmesh, par, 'exp_graded', 1);
dev.gradNv = build_property(par.Nv, xmesh, par, 'exp_graded', 1);

% Surface recombination velocity equivalence schemes
dev.taun_vsr = build_property(par.taun, xmesh, par, 'taun_vsr', 0);
dev.taup_vsr = build_property(par.taup, xmesh, par, 'taup_vsr', 0);
dev.alpha0 = build_property(0, xmesh, par, 'alpha0', 1);
dev.beta0 = build_property(0, xmesh, par, 'beta0', 1);
dev.alpha0_xn = build_property(0, xmesh, par, 'alpha0_xn', 1);
dev.beta0_xp = build_property(0, xmesh, par, 'beta0_xp', 1);
dev.dint = build_property(0, xmesh, par, 'dint', 1);

% Switches
dev.int_switch = build_property(0, xmesh, par, 'int_switch', 1);
dev.bulk_switch = abs(dev.int_switch-1);
if par.vsr_mode
    dev.vsr_zone = build_property(0, xmesh, par, 'vsr_zone', 1);
    dev.srh_zone = dev.bulk_switch;
    dev.Field_switch = dev.bulk_switch;
    
    dev.NA = build_property(par.NA, xmesh, par, 'zeroed', 0);
    dev.ND = build_property(par.ND, xmesh, par, 'zeroed', 0);
    dev.epp = build_property(par.epp, xmesh, par, 'constant', 0); 
    dev.ni = build_property(par.ni, xmesh, par, 'constant', 0);
    dev.nt = build_property(par.nt, xmesh, par, 'constant', 0);
    dev.pt = build_property(par.pt, xmesh, par, 'constant', 0);
else 
    dev.vsr_zone = zeros(1, length(xmesh));
    dev.srh_zone = ones(1, length(xmesh));
    dev.Field_switch = ones(1, length(xmesh));
    
    dev.NA = build_property(par.NA, xmesh, par, 'exp_graded', 0);
    dev.ND = build_property(par.ND, xmesh, par, 'exp_graded', 0);
    dev.epp = build_property(par.epp, xmesh, par, 'lin_graded', 0);
    dev.ni = build_property(par.ni, xmesh, par, 'exp_graded', 0);
    dev.nt = build_property(par.nt, xmesh, par, 'exp_graded', 0);
    dev.pt = build_property(par.pt, xmesh, par, 'exp_graded', 0);
end

% Tranlsated co-ordinates
dev.xprime = build_property(par.xx, xmesh, par, 'xprime', 1);
dev.xprime_n = build_property(par.xx, xmesh, par, 'xprime_n', 1);
dev.xprime_p = build_property(par.xx, xmesh, par, 'xprime_p', 1);
dev.sign_xn = build_property(par.xx, xmesh, par, 'sign_xn', 1);
dev.sign_xp = build_property(par.xx, xmesh, par, 'sign_xp', 1);

end
