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
dcum = par.dcum;

switch par.xmesh_type
    % Linearly spaced
    case 1
        x = linspace(0,dcum(end),par.parr(1)+par.parr(2)+par.parr(3));
    case 2
        %x = zeros(1, 600);
        parrcum = [0, cumsum(par.parr)];
        parr = par.parr;
        dcum = [0, par.dcum];
        j = 1;
        k = 1;
        % Iteration around number of layers and desired points
        for i=1:2*length(par.parr)-1
            
            if rem(i, 2) == 1
                parrint(j) = parr(k);
                j = j+1;
                k = k+1;
            elseif rem(i, 2) == 0
                parrint(j) = par.pint;
                j = j+1;
            end
        end
        parrintcum = cumsum(parrint);
        parrintcum = [0,parrintcum];
        j=1;
        k=1;
        
        for i=1:length(parrint)
            if i == 1
                A = 0;
            else
                A = 1;
            end
            
            if i == length(parrint)
                B = 0;
            else
                B = 1;
            end
            
            if rem(i, 2) == 1
                linarr = linspace(dcum(j)+A*(par.dint+par.dint/par.pint), dcum(j+1)-B*(par.dint+par.dint/par.pint), parrint(i));
                x(1, (parrintcum(i)+1):parrintcum(i+1)) = linarr;
            elseif rem(i, 2) == 0
                linarr = linspace(dcum(j+1)-par.dint, dcum(j+1)+par.dint, parrint(i));
                x(1, (parrintcum(i)+1):parrintcum(i+1)) = linarr;
                j=j+1;
            end
            
            
        end
        
        
    case 3
        % build number array from dcell (cell array containing lengths)
        dcellarr = 0;
        for i = 1:size(par.dcell, 1)
            p_temp = par.dcell{i,:};
            dcellarr = [dcellarr, cell2mat(p_temp)];
        end
        dcellcum = cumsum(dcellarr);
        
        % build number array from layer_points (cell array containing lengths)
        p_array = 0;
        for i = 1:size(par.layer_points, 1)
            p_temp = par.layer_points{i,:};
            p_array = [p_arr, cell2mat(p_temp)];
        end
        p_cum = cumsum(p_array);
        
        j = 1;
        k = 1;
        % Iteration around number of layers and desired points
        for i=1:2*length(par.d)-1
            % i tracks the stack layers including interfaces
            % j tracks the point layer (not inc interfaces)
            % k tracks the stack layers not including interfaces
            
            if rem(i, 2) == 1
                for n=1:length(par.layer_points{k,:})
                    p_temp = par.layer_points{k, :};
                    stacklayerarr = cell2mat(p_temp);
                    parrint(j) = stacklayerarr(n);
                    
                    p_temp = par.dcell{k, :};
                    stacklayerarr = cell2mat(p_temp);
                    darrint(j) = stacklayerarr(n);
                    
                    j = j+1;
                end
                k = k+1;
            elseif rem(i, 2) == 0
                parrint(j) = par.pint;
                darrint(j) = par.dint;
                j = j+1;
            end
            
        end
        parrintcum = cumsum(parrint);
        parrintcum = [0,parrintcum];
        darrintcum = cumsum(darrint);
        darrintcum = [0,darrintcum];
        
        for i=1:length(parrint)     % Layer with interface index
            if i == 1
                A = 0;
            else
                A = 1;
            end
            
            if i == length(parrint)
                B = 0;
            else
                B = 1;
            end
            
            % Write the entries into x
            linarr = linspace(darrintcum(i), darrintcum(i+1)-B*(darrint(i)/parrint(i)), parrint(i));
            x(1, (parrintcum(i)+1):parrintcum(i+1)) = linarr;
            
        end
        % To account for rounding errors
        x(end) = darrintcum(end);
        
    case 4
        d = par.dcell;
        p = par.layer_points;
        dcum0 = par.dcum0;
        
        xcell = cell(1,length(p));
        for i=1:length(par.layer_points)
            linarr = linspace(dcum0(i), dcum0(i+1)-(d(i)/p(i)), p(i));
            xcell{i} = linarr;
        end
        x = [xcell{:}];
        x = [x, dcum0(end)];
        
    case 5
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
        for i = 1:length(p)
            if any(strcmp(par.layer_type{1,i}, {'layer', 'active'}))
                parr = -0.5 : (1/p(i)) : 0.5;
                x_layer = erf(2*pi*xmesh_coeff(i)*parr);
                x_layer = normalize(x_layer,'range');
                x_layer = dcum0(i) + x_layer * d(i);
                xcell{i} = x_layer(1:end-1);
            elseif strcmp(par.layer_type{1,i}, 'junction')
                x_layer = linspace(dcum0(i), dcum0(i+1)-(d(i)/p(i)), p(i));
                xcell{i} = x_layer;
            end
        end
        x = [xcell{:}];
        x = [x, dcum0(end)];

    otherwise
        error('DrIFtFUSION:xmesh_type', [mfilename ' - xmesh_type not recognized'])
end

px = length(x);

if par.meshx_figon == 1
    
    parr = 1:1:px;
    dpdx = gradient(parr, x);
    xmir = x;
    pxmir = 1:1:length(x);
    
    figure(1010);
    plot(xmir, pxmir,'.')
    xlabel('Position');
    ylabel('Point');
    
    figure(1011)
    plot(xmir, dpdx);
    xlabel('Position');
    ylabel('dpdx');
    
    drawnow
end

end