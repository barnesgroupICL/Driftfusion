function [x] = meshgen_x(params)

v2struct(params);

meshfigon = 1; % concentration of mesh points image

% Linearly spaced
if xmesh_type == 1

    x = linspace(0,xmax,pp+pii+pn);

% Linearly spaced, more points at interfaces   
elseif xmesh_type == 2
   
    x = [linspace(0, tp-tinter, pp),...
        linspace(tp-tinter+deltax, tp, pinter),...
        linspace(tp+deltax, tp+tinter, pinter),...
        linspace(tp+tinter+deltax, tp+ti-tinter-deltax, pii),...
        linspace(tp+ti-tinter, tp+ti-deltax, pinter),...
        linspace(tp+ti, tp+ti+tinter-deltax, pinter),...
        linspace(tp+ti+tinter, xmax, pn)]; 
    
% Linearly spaced, more points at interfaces and electrodes
elseif xmesh_type == 3
   
    x = [linspace(0, te, pepe),...
        linspace(te+deltax, tp-tinter, pp),...
        linspace(tp-tinter+deltax, tp, pinter),...
        linspace(tp+deltax, tp+tinter, pinter),...
        linspace(tp+tinter+deltax, tp+ti-tinter-deltax, pii),...
        linspace(tp+ti-tinter, tp+ti-deltax, pinter),...
        linspace(tp+ti, tp+ti+tinter-deltax, pinter),...
        linspace(tp+ti+tinter, xmax-te-deltax, pn),...
        linspace(xmax-te, xmax, pepe);]; 
        
end 

px = length(x);

if meshfigon == 1

    xmir = x;
    pxmir = 1:1:length(x);
    
    figure(1010);
    plot(xmir, pxmir, '.');
    xlabel('Position');
    ylabel('Point');

end

end
