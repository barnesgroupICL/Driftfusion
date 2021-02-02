function [t] = meshgen_t(par)
% Creates the output time mesh
% Note- ODE15S uses an adaptive time step and then interpolates to the
% requested times
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
switch par.tmesh_type
    case 1
        t = linspace(0,par.tmax,par.tpoints);     
    case 2
        t = logspace(log10(par.t0),log10(par.tmax),par.tpoints);
        %To avoid rounding errors
        t(1) = 0;
        t(end) = par.tmax;
    case 3
        % 2 log meshes consecutively
        t1 = logspace(log10(par.t0),log10(par.tmax/2),round(par.tpoints/2));
        t1(1) = 0;
        t1(end) = par.tmax/2;
        t2 = t1(end) + logspace(log10(par.t0),log10(par.tmax/2),round(par.tpoints/2));
        t2(end) = par.tmax;
        t = [t1, t2];
end

if par.mesht_figon == 1
    
    ptmir = 1:1:length(t);
    
    figure(200);
    plot(t, ptmir, '.');
    xlabel('Position');
    ylabel('Time'); 
end

end