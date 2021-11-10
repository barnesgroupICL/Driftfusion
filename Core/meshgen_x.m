function x = meshgen_x(par)
% Generates the spatial mesh dependent on option defined by XMESH_TYPE
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
xmesh_type = par.xmesh_type;

switch xmesh_type
    % Linearly spaced
    case 'linear'
        d = par.dcell;
        p = par.layer_points;
        dcum0 = par.dcum0;
        
        xcell = cell(1,length(p));
        for i = 1:length(par.layer_points)
            linarr = linspace(dcum0(i), dcum0(i+1)-(d(i)/p(i)), p(i));
            xcell{i} = linarr;
        end
        x = [xcell{:}];
        x = [x, dcum0(end)];
        
    case 'erf-linear'
        % Error function for each layer
        d = par.dcell;
        p = par.layer_points;
        dcum0 = par.dcum0;
        
        % For backwards compatibility
        if length(par.xmesh_coeff) < length(p)
            xmesh_coeff = 0.7*ones(1, length(p));
        else
            xmesh_coeff = par.xmesh_coeff;
        end
        
        xcell = cell(1,length(p)); 
        for i = 1:length(par.layer_type)
            if any(strcmp(par.layer_type{1,i}, {'layer', 'active'}))
                parr = -0.5 : (1/p(i)) : 0.5;
                x_layer = erf(2*pi*xmesh_coeff(i)*parr);
                x_layer = x_layer - x_layer(1);       % Subtract base to get zero
                x_layer  = x_layer./max(x_layer);   % Normalise the funciton
                x_layer = dcum0(i) + x_layer * d(i);
                xcell{i} = x_layer(1:end-1);
            elseif any(strcmp(par.layer_type{1,i}, {'junction', 'interface'}))
                x_layer = linspace(dcum0(i), dcum0(i+1)-(d(i)/p(i)), p(i));
                xcell{i} = x_layer;
            end
        end
        x = [xcell{:}];
        x = [x, dcum0(end)];

    otherwise
        error('xmesh_type not recognized')
end

px = length(x);
end