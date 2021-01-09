function devprop = build_property(property, xmesh, par, interface_switch, gradient_property)
% Builds the device property - i.e. defines the properties at
% every x position in the device
% PROPERTY          - the variable name of the propery e.g. par.EA
% XMESH             - as name suggests
% PAR               - parameters object
% INTERFACE_SWICTH  -   'zeroed' = set property value to zero for interfaces
%                       'constant' = constant property values in interfaces
%                       'lin_graded' = graded property values in interfaces
%                       'log_graded' = graded property values in interfaces
% GRADIENT_PROPERTIES - 1 if the property is a gradient e.g. dEAdx, 0 otherwise
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
devprop = zeros(1, length(xmesh));

for i=1:length(par.dcum)                % i is the layer index
    for j = 1:length(xmesh)             % j is the spatial mesh point index
        if any(strcmp(par.layer_type{1,i}, {'layer', 'active'})) == 1
            if xmesh(j) >= par.dcum0(i)
                if gradient_property == 1
                    devprop(j) = 0;
                else
                    devprop(j) = property(i);
                end
            end
        elseif any(strcmp(par.layer_type{1,i}, {'junction'})) == 1
            if xmesh(j) >= par.dcum0(i)
                xprime = xmesh(j)-par.dcum0(i);
                if xmesh(j) >= par.dcum0(i)
                    deff = par.d(i);
                    switch interface_switch
                        case 'zeroed'
                            devprop(j) = 0;
                        case 'constant'
                            devprop(j) = property(i);
                        case 'lin_graded'
                            gradient = (property(i+1)-property(i-1))/deff;
                            if gradient_property == 1
                                devprop(j) = gradient;
                            else
                                devprop(j) = property(i-1) + xprime*gradient;
                            end
                        case 'log_graded'
                            log_gradient = (log(property(i+1))-log(property(i-1)))/deff;
                            if gradient_property == 1
                                devprop(j) = property(i-1)*log_gradient*exp(log_gradient*xprime);
                            else
                                devprop(j) = property(i-1)*exp(log_gradient*xprime);
                            end
                        case 'surface_rec_taun'
                            alpha = (par.EA(i-1) - par.EA(i+1))/(par.kB*par.T) + log(par.Nc(i+1))-log(par.Nc(i-1));
                            if alpha <= 0
                                devprop(j) = (deff/par.sn(i))*exp(alpha*xprime/deff);
                            elseif alpha > 0
                                devprop(j) = (deff/par.sn(i))*exp(alpha*(xprime-deff)/deff);
                            end
                        case 'surface_rec_taup'
                            beta = (par.IP(i+1) - par.IP(i-1))/(par.kB*par.T) + log(par.Nv(i+1))-log(par.Nv(i-1));
                            if beta < 0
                                devprop(j) = (deff/par.sp(i))*exp(beta*xprime/deff);
                            elseif beta > 0
                                devprop(j) = (deff/par.sp(i))*exp(beta*(xprime-deff)/deff);
                            end
                        case 'surface_rec_nt'
                            alpha = (par.EA(i-1) - par.EA(i+1))/(par.kB*par.T) + log(par.Nc(i+1))-log(par.Nc(i-1));
                            if alpha <= 0
                                 devprop(j) = par.nt(i)*exp(alpha*xprime/deff);
                            elseif alpha > 0
                                 devprop(j) = par.nt(i)*exp(alpha*(xprime-deff)/deff);
                            end      
                        case 'surface_rec_pt'
                            beta = (par.IP(i+1) - par.IP(i-1))/(par.kB*par.T) + log(par.Nv(i+1))-log(par.Nv(i-1)); 
                            if beta < 0
                                devprop(j) = par.pt(i)*exp(beta*xprime/deff);
                            elseif beta > 0
                                devprop(j) = par.pt(i)*exp(beta*(xprime-deff)/deff);
                            end 
                    end
                end
            end
        end
    end
end
end
