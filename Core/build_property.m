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
                    end
                end
            end
        end
    end
end
end
