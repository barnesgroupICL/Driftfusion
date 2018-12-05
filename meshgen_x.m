function x = meshgen_x(par)

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
        
%         parrcum = [0, cumsum(par.parr)];
%         parr = par.parr;
%         darr = par.d;
%         dcum = [0, par.dcum];
        
        % build number array from dcell (cell array containing lengths)
        dcellarr = 0;
        for i = 1:size(par.dcell, 1)
            tempcell = par.dcell{i,:};
            dcellarr = [dcellarr, cell2mat(tempcell)];
        end
        dcellcum = cumsum(dcellarr);
        
        % build number array from pcell (cell array containing lengths)
        pcellarr = 0;
        for i = 1:size(par.pcell, 1)
            tempcell = par.pcell{i,:};
            pcellarr = [pcellarr, cell2mat(tempcell)];
        end
        pcellcum = cumsum(pcellarr);
        
        j = 1;
        k = 1;
        % Iteration around number of layers and desired points
        for i=1:2*length(par.d)-1
            % i tracks the stack layers including interfaces
            % j tracks the point layer (not inc interfaces)
            % k tracks the stack layers not including interfaces
            
            if rem(i, 2) == 1
                for n=1:length(par.pcell{k,:})
                    tempcell = par.pcell{k, :};
                    stacklayerarr = cell2mat(tempcell);
                    parrint(j) = stacklayerarr(n);
                    
                    tempcell = par.dcell{k, :};
                    stacklayerarr = cell2mat(tempcell);
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
    %ylim([-10,10])
    
    drawnow
end

end