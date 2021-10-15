function devprop = build_property(property, xmesh, par, interface_switch, gradient_property)
% Builds the device property - i.e. defines the properties at
% every x position in the device
% PROPERTY          - the variable name of the propery e.g. par.Phi_EA
% XMESH             - as name suggests
% PAR               - parameters object
% INTERFACE_SWICTH  -   'zeroed' = set property value to zero for interfaces
%                       'constant' = constant property values in interfaces
%                       'lin_graded' = graded property values in interfaces
%                       'exp_graded' = graded property values in interfaces
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
                % Gradient coefficient for surface recombination equivalence
                alpha0 = ((par.Phi_EA(i-1) - par.Phi_EA(i+1))/(par.kB*par.T) + (log(par.Nc(i+1))-log(par.Nc(i-1))))/deff;
                if alpha0 < 0
                    xprime_n = xprime;
                    alpha0_xn = alpha0;     % the sign of alpha0 is referenced to the direction of xprime_n
                else
                    xprime_n = deff-xprime;
                    alpha0_xn = -alpha0;     % the sign of alpha0 is referenced to the direction of xprime_n
                end
                % Gradient coefficient for surface recombination equivalence
                beta0 = ((par.Phi_IP(i+1) - par.Phi_IP(i-1))/(par.kB*par.T) + (log(par.Nv(i+1))-log(par.Nv(i-1))))/deff;
                if beta0 < 0
                    xprime_p = xprime;
                    beta0_xp = beta0;        % the sign of beta is referenced to the direction of xprime_p
                else
                    xprime_p = deff-xprime;  
                    beta0_xp = -beta0;       % the sign of beta is referenced to the direction of xprime_p
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
                    case 'exp_graded'
                        log_gradient = (log(property(i+1))-log(property(i-1)))/deff;
                        if gradient_property == 1
                            devprop(j) = property(i-1)*log_gradient*exp(log_gradient*xprime);
                        else
                            devprop(j) = property(i-1)*exp(log_gradient*xprime);
                        end
                    case 'taun_vsr'
                        if par.sn(i) == 0
                            devprop(j) = 1e100; % Avoid divide by zero error
                        else
                            devprop(j) = deff*par.frac_vsr_zone/(par.sn(i));
                        end
                    case 'taup_vsr'
                        if par.sp(i) == 0
                            devprop(j) = 1e100; % Avoid divide by zero error
                        else
                            devprop(j) = deff*par.frac_vsr_zone/(par.sp(i));
                        end
                    case 'mu_n_vsr'
                        if alpha0 < 0
                            devprop(j) = par.mu_n(i-1)*exp(abs(alpha0)*xprime_n);
                        elseif alpha0 > 0
                            devprop(j) = par.mu_n(i+1)*exp(abs(alpha0)*xprime_n);
                        else
                            devprop(j) = par.mu_n(i);
                        end
                    case 'mu_p_vsr'
                        if beta0 < 0
                            devprop(j) = par.mu_p(i-1)*exp(abs(beta0)*xprime_p);
                        elseif beta0 > 0
                            devprop(j) = par.mu_p(i+1)*exp(abs(beta0)*xprime_p);
                        else
                            devprop(j) = par.mu_p(i);
                        end
                    case 'int_switch'
                        devprop(j) = 1;
                    case 'xprime'
                        devprop(j) = xprime;
                    case 'xprime_n'
                        devprop(j) = xprime_n; 
                    case 'xprime_p'
                        devprop(j) = xprime_p;
                    case 'sign_xn'
                        if alpha0 < 0
                            devprop(j) = 1;
                        else
                            devprop(j) = -1;
                        end
                    case 'sign_xp'
                        if beta0 < 0
                            devprop(j) = 1;
                        else
                            devprop(j) = -1;
                        end
                    case 'dint'
                        devprop(j) = deff;
                    case 'alpha0'
                        devprop(j) = alpha0;
                    case 'beta0'
                        devprop(j) = beta0;
                    case 'alpha0_xn'
                        devprop(j) = alpha0_xn;
                    case 'beta0_xp'
                        devprop(j) = beta0_xp;
                    case 'vsr_zone'
                        if par.vsr_zone_loc(i) == "L"
                            if xprime <= deff*par.frac_vsr_zone
                                devprop(j) = 1;
                            else
                                devprop(j) = 0;
                            end
                        elseif par.vsr_zone_loc(i) == "R"
                            if xprime >= deff*(1 - par.frac_vsr_zone) && xprime <= deff
                                devprop(j) = 1;
                            else
                                devprop(j) = 0;
                            end
                        elseif par.vsr_zone_loc(i) == "C"
                            if xprime >= (deff/2)*(1 - par.frac_vsr_zone) && xprime <= (deff/2)*(1 + par.frac_vsr_zone)
                                devprop(j) = 1;
                            else
                                devprop(j) = 0;
                            end
                        end
                end
            end
        end
    end
end
end

