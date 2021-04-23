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
            if i == 1
                condition = (xmesh(j) >= par.dcum0(i));
            else
                condition = (xmesh(j) >= par.dcum0(i));
            end
            
            if condition
                if gradient_property == 1
                    devprop(j) = 0;
                else
                    devprop(j) = property(i);
                end
            end
        elseif any(strcmp(par.layer_type{1,i}, {'junction', 'interface'})) == 1

            if xmesh(j) >= par.dcum0(i)
                xprime = xmesh(j)-par.dcum0(i);
                    deff = par.d(i);
                    % Gradient coefficients for surface recombination equivalence
                    alpha = ((par.EA(i-1) - par.EA(i+1))/(par.kB*par.T) + (log(par.Nc(i+1))-log(par.Nc(i-1))))/deff;
                    if alpha < 0
                        xprime_n = xprime;
                    else
                        xprime_n = deff-xprime;
                    end

                    beta = ((par.IP(i+1) - par.IP(i-1))/(par.kB*par.T) + (log(par.Nv(i+1))-log(par.Nv(i-1))))/deff;
                    if beta < 0
                        xprime_p = xprime;
                    else
                        xprime_p = deff-xprime;
                    end

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
                        case 'taun_vsr'
                                devprop(j) = (deff/par.sn(i));%*exp(-abs(alpha)*xprime_n);
                        case 'taup_vsr'
                                devprop(j) = (deff/par.sp(i));%*exp(-abs(beta)*xprime_p);
                        case 'nt_vsr'
                                devprop(j) = par.nt(i);%*exp(-abs(alpha)*xprime_n);
                        case 'pt_vsr'
                                devprop(j) = par.pt(i);%*exp(-abs(beta)*xprime_p);
                        case 'ni_vsr'
                                devprop(j) = par.ni(i);%./abs((exp(abs(alpha)*xprime_n).*exp(abs(beta)*xprime_p)).^0.5);
                        case 'mue_interface'
                            if alpha < 0
                                devprop(j) = (par.mue(i-1)*exp(abs(alpha)*xprime_n));%(par.mue(i-1)+((par.Nv(i+1)/par.Nv(i+1))*(par.sn(i) + par.sp(i))*(1/(par.kB*par.T*abs(alpha)))))*exp(abs(alpha)*xprime_n);
                                %devprop(j) = par.mue(i-1)+(((par.Nv(i+1)/par.Nv(i+1))*(par.sn(i) + par.sp(i))*(1/(par.kB*par.T*abs(alpha)))))*exp(abs(alpha)*xprime_n);
                            elseif alpha > 0
                                devprop(j) = (par.mue(i+1)*exp(abs(alpha)*xprime_n));%(par.mue(i+1)+((par.Nv(i-1)/par.Nv(i-1))*(par.sn(i) + par.sp(i))*(1/(par.kB*par.T*abs(alpha)))))*exp(abs(alpha)*xprime_n);
                                %devprop(j) = par.mue(i+1)+(((par.Nv(i-1)/par.Nv(i-1))*(par.sn(i) + par.sp(i))*(1/(par.kB*par.T*abs(alpha)))))*exp(abs(alpha)*xprime_n);
                            else
                                devprop(j) = max([par.mue(i-1),par.mue(i+1)]);
                            end
                        case 'muh_interface'
                            if beta < 0
                                devprop(j) = (par.muh(i-1)*exp(abs(beta)*xprime_p));%(par.muh(i-1)+((par.Nc(i+1)/par.Nv(i-1))*(par.sn(i) + par.sp(i))*(1/(par.kB*par.T*abs(beta)))))*exp(abs(beta)*xprime_p);
                                %devprop(j) = par.muh(i-1)+(((par.Nc(i+1)/par.Nv(i-1))*(par.sn(i) + par.sp(i))*(1/(par.kB*par.T*abs(beta)))))*exp(abs(beta)*xprime_p);
                            elseif beta > 0
                                devprop(j) = (par.muh(i+1)*exp(abs(beta)*xprime_p));%(par.muh(i+1)+((par.Nc(i-1)/par.Nv(i+1))*(par.sn(i) + par.sp(i))*(1/(par.kB*par.T*abs(beta)))))*exp(abs(beta)*xprime_p);
                                %devprop(j) = par.muh(i+1)+(((par.Nc(i-1)/par.Nv(i+1))*(par.sn(i) + par.sp(i))*(1/(par.kB*par.T*abs(beta)))))*exp(abs(beta)*xprime_p);
                            else
                                devprop(j) = max([par.muh(i-1),par.muh(i+1)]);
                            end
                        case 'int_switch'
                            devprop(j) = 1;
                        case 'xprime_n'
                            devprop(j) = xprime_n;
                        case 'xprime_p'
                            devprop(j) = xprime_p;
                        case 'xprime'
                            devprop(j) = xprime_p;
                        case 'alpha'
                            devprop(j) = alpha;
                        case 'beta'
                            devprop(j) = beta;
                    end
            end
        end
    end
end
end
