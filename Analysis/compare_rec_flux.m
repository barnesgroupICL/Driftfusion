function sigma_sum_filter = compare_rec_flux(sol, RelTol_vsr, AbsTol_vsr, plot_switch)
% Script to compare the interfacial recombination fluxes from Driftfusion.
% The integrated recombination flux from the volumetric surface
% recombination model is compared with the flux calculated from the
% 2D SRH expression using the boundary carrier densities, ns and ps
% directly
%
% The sum of the fractional differences SIGMA_SUM is calculated by first
% summing the recombination fluxes for all interfaces and then calculating
% the fraction.
%
%% LICENSE
% Copyright (C) 2021  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Input arguments
% PLOT_SWITCH = 1 for plot outputs
% RELTOL_VSR is the fractional difference between DF 2D and volumetric calculations above
% which a warning is displayed
% ABSTOL_VSR is the absolute value of the recombination flux below which
% the difference is ignored (lower values tend to have much higher relative
% errors)
%
%% Start code
par = sol.par;
u = sol.u;
t = sol.t;
pcum1 = par.pcum + 1;   % Cumulative layer points array
dev = par.dev_sub;
x = par.xx;
x_sub = par.x_sub;

% Calculate alpha and beta- note these are referenced to the direction of x
% NOT x_n and x_p
dVdx = zeros(length(t), length(x_sub));
for i = 1:length(t)
    dVdx(i,:) = (u(i, 2:end, 1) - u(i, 1:end-1, 1))./(x(2:end) - x(1:end-1));
end
alpha0 = repmat(dev.alpha0, length(t), 1);
beta0 = repmat(dev.beta0, length(t), 1);
alpha = par.q.*dVdx./(par.kB*par.T) + alpha0;
beta = par.q.*-dVdx./(par.kB*par.T) + beta0;

n = u(:, :, 2);
p = u(:, :, 3);

%% Driftfusion 3D
rx = dfana.calcr(sol, "sub");

%% Get indexes of interfaces
%int_index = find(contains(par.layer_type, 'interface'));   % only
%compatible from 2016 onwards
int_logical = zeros(1, length(par.layer_type));
for i = 1:length(par.layer_type)
    int_logical(i) = any(any(strcmp(par.layer_type(i), {'interface', 'junction'})));
end
loc = find(int_logical); % interface layer locations

ns = zeros(length(t), length(loc));     % Store time array of each ns in new column
ps = zeros(length(t), length(loc));     % Store time array of each ps in new column
R_abrupt = zeros(length(t), length(loc));
R_vsr = zeros(length(t), length(loc));
sigma = zeros(length(t), length(loc));
legstr_R = cell(2*length(loc) + 1, 1);
legstr_sigma = cell(length(loc) + 1, 1);

for i = 1:length(loc)
    sn = par.sn(loc(i));
    sp = par.sp(loc(i));
    ni = par.ni(loc(i));
    nt = par.nt(loc(i));
    pt = par.pt(loc(i));
    
    % Check location of interfacial surface carrier densities
    p_L = pcum1(loc(i)-1);
    p_R = pcum1(loc(i));
    
    for k = 1:length(t)
        if alpha(k, p_L) <= 0
            ns(k, i) = n(k, p_L);
        elseif alpha(k, p_L) > 0
            ns(k, i) = n(k, p_R);
        end
        
        if beta(k, p_L) <= 0
            ps(k, i) = p(k, p_L);
        elseif beta(k, p_L) > 0
            ps(k, i) = p(k, p_R);
        end
        R_vsr(k, i) = trapz(x_sub(p_L-1:p_R+1), rx.vsr(k, p_L-1:p_R+1), 2);
    end
    
    R_abrupt(:, i) = (ns(:, i).*ps(:, i) - ni^2)./((1/sn).*(ps(:, i) + pt) + (1/sp).*(ns(:, i) + nt));
    
    %% Fractional difference
    sigma(:, i) = ((R_abrupt(:, i) - R_vsr(:, i))./R_abrupt(:, i));
        
    if plot_switch
        figure(300)
        semilogy(t, R_abrupt(:, i) , t, R_vsr(:, i) , '-.')
        xlabel('Time [s]')
        ylabel('Recombination flux [cm-2s-1]')
        xlim([t(1), t(end)])
        legstr_R{2*i - 1} = ['Interface', num2str(i), '- abrupt'];
        legstr_R{2*i} = ['Interface', num2str(i), '- discrete'];
        hold on
        
        figure(301)
        semilogy(t, ns(:, i) , t, ps(:, i))
        xlabel('Time [s]')
        ylabel('Carrier density')
        legstr_nsps{2*i - 1} = ['Interface', num2str(i), '- ns'];
        legstr_nsps{2*i} = ['Interface', num2str(i), '- ps'];
        xlim([t(1), t(end)])
        hold on
        
        figure(302)
        plot(t, sigma(:, i) )
        xlabel('Time [s]')
        ylabel('Fractional difference')
        xlim([t(1), t(end)])
        legstr_sigma{i} = ['Interface', num2str(i)];
        hold on
    end
end

R_vsr_sum = sum(R_vsr, 2);
R_vsr_filter = R_vsr_sum;
R_vsr_filter(R_vsr_filter < AbsTol_vsr) = NaN;
R_abrupt_sum = sum(R_abrupt, 2);
sigma_sum = (1-(R_vsr_sum./R_abrupt_sum));
sigma_sum_filter = (1-(R_vsr_filter./R_abrupt_sum));

if max(abs(sigma_sum_filter)) > RelTol_vsr
    warning(['The max volumetric surface recombination model fractional error (sigma_max = ', num2str(max(abs(sigma_sum_filter))),') for recombination fluxes above ', num2str(AbsTol_vsr), ' cm-2s-1 exceeded the user-defined tolerance level (tol_vsr = ', ...
        num2str(RelTol_vsr), '). Consider:' newline...
        '- increasing the interface layer electronic mobilities' newline...
        '- reducing energetic barriers' newline...
        '- manually relocating recombination zones' newline...
        '- reducing interfacial surface recombination velocities' newline...
        '- switching to an alternative recombination model'])
end
    
if plot_switch
    figure(300)
    semilogy(t, AbsTol_vsr*ones(1, length(t)), 'k--')
    legstr_R{end} = 'AbsTol';
    legend(legstr_R);
    
    figure(301)
    legend(legstr_nsps)
    
    figure(302)
    plot(t, sigma_sum_filter, 'k--')
    xlabel('Time [s]')
    ylabel('Absolute fractional difference')
    xlim([t(1), t(end)])
    legstr_sigma{end} = 'Sum (filtered)';
    legend(legstr_sigma)
    
    figure(300); hold off
    figure(301); hold off
    figure(302); hold off
end
end