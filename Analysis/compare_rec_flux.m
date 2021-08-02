function compare_rec_flux(sol_df, RelTol_vsr, AbsTol_vsr, plot_switch)
% Script to compare the interfacial recombination fluxes from Driftfusion
% (DF) and IonMonger (IM)
% Currently only working for ETL-AL-HTL architecture
% PLOT_SWITCH = 1 for plot outputs
% RELTOL_VSR is the fractional difference between DF 2D and volumetric calculations above
% which a warning is displayed
% ABSTOL_VSR is the absolute value of the recombination flux below which
% the difference is ignored (lower values tend to have much higher relative
% errors)

par = sol_df.par;
u = sol_df.u;
t = sol_df.t;
x = sol_df.x;
al = par.active_layer;
pcum1 = par.pcum + 1;   % Cumulative layer points array
dev = par.dev;
alpha0 = dev.alpha0;
beta0 = dev.beta0;
x = par.x;
xsub = par.x_ihalf;

n = u(:, :, 2);
p = u(:, :, 3);

%% Driftfusion 3D
rx = dfana.calcr_ihalf(sol_df);

%% Get indexes of interfaces
%int_index = find(contains(par.layer_type, 'interface'));   % only
%compatible from 2016 onwards
int_logical = strcmp(par.layer_type, 'interface');
loc = find(int_logical); % interface layer locations

ns = zeros(length(t), length(loc));     % Store time array of each ns in new column
ps = zeros(length(t), length(loc));     % Store time array of each ps in new column
R_abrupt = zeros(length(t), length(loc));
R_vsr = zeros(length(t), length(loc));
R_vsr_filter = NaN(length(t), length(loc));
R_abrupt_filter = NaN(length(t), length(loc));
sigma = zeros(length(t), length(loc));

for i = 1:length(loc)
    sn = par.sn(loc(i));
    sp = par.sp(loc(i));
    ni = par.ni(loc(i));
    nt = par.nt(loc(i));
    pt = par.pt(loc(i));
    
    % Check location of interfacial surface carrier densities
    if alpha0(pcum1(loc(i)-1)) <= 0
        ns(:, i) = n(:, pcum1(loc(i)-1));
    elseif alpha0(pcum1(loc(i)-1)) > 0
        ns(:, i) = n(:, pcum1(loc(i)));
    end
    
    if beta0(pcum1(loc(i)-1)) <= 0
        ps(:, i) = p(:, pcum1(loc(i)-1));
    elseif beta0(pcum1(loc(i)-1)) > 0
        ps(:, i) = p(:, pcum1(loc(i)));
    end
    
    R_abrupt(:, i) = (ns(:, i).*ps(:, i) - ni^2)./((1/sn).*(ps(:, i) + pt) + (1/sp).*(ns(:, i) + nt));
    
    for k = 1:length(t)
        R_vsr(k, i) = trapz(xsub(pcum1(loc(i)-1)-1:pcum1(loc(i))), rx.vsr(k, pcum1(loc(i)-1)-1:pcum1(loc(i))), 2);
    end
    
    %% Fractional difference
    sigma(:, i) = (R_abrupt(:, i)-R_vsr(:, i))./R_abrupt(:, i);
    
    R_vsr_filter(find(R_vsr(:, i) < AbsTol_vsr), i) = NaN;
    R_abrupt_filter(find(R_abrupt(:, i) < AbsTol_vsr), i) = NaN;
        
    if plot_switch
        figure(300)
        semilogy(t, R_abrupt(:, i) , t, R_vsr(:, i) , 'k--')
        xlabel('Time [s]')
        ylabel('Recombination flux [cm-2s-1]')
        legend('Abrupt interface', 'Discrete interface (volumetric)')
        xlim([t(1), t(end)])
        hold on
        
        figure(301)
        semilogy(t, ns(:, i) , t, ps(:, i))
        xlabel('Time [s]')
        ylabel('Carrier density')
        legend('ns', 'ps')
        xlim([t(1), t(end)])
        hold on
        
        figure(302)
        plot(t, sigma(:, i) )
        xlabel('Time [s]')
        ylabel('Fractional difference')
        xlim([t(1), t(end)])
        hold on
    end
end

sigma_sum = 1-(sum(R_vsr, 2)./sum(R_abrupt, 2));
sigma_sum_filter = sigma_sum;
[AbsTol_row, AbsTol_col] = find(R_vsr < AbsTol_vsr);
sigma_sum_filter(AbsTol_row) = NaN;

if max(abs(sigma_sum_filter)) > RelTol_vsr
    warning(['The max volumetric surface recombination model fractional error (sigma_max = ', num2str(max(abs(sigma_sum))),') for recombination fluxes above ', num2str(AbsTol_vsr), ' cm-2s-1 exceeded the user-defined tolerance level (tol_vsr = ', ...
        num2str(RelTol_vsr), '). Consider:' newline...
        '- increasing the interface layer electronic mobilities' newline...
        '- reducing energetic barriers' newline...
        '- manually relocating recombination zones' newline...
        '- reducing interfacial surface recombination velocities' newline...
        '- switching to an alternative recombination model'])
end
    
if plot_switch
    figure(303)
    plot(t, sigma_sum_filter)
    xlabel('Time [s]')
    ylabel('Fractional difference')
    figure(300); hold off
    figure(301); hold off
    figure(302); hold off
end

end
