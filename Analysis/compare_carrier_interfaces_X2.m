function compare_carrier_interfaces_X2(sol, tarr)
% Function to compare the interfacial carrier densities between simulation
% and analytical solution
% v5 Calculated on the whole mesh, uses jn instead of js
% substitutes r for r=(js-j)/x
% v6 Tests expressions for ns, ps for majority carrier sides of the
% v7 as with v6 but calculated on the half mesh
% v8 uses left hand boundary values only and adapts expression for n(d) =
% ns and p(d) = ps dependent on sign of alpha and beta
% v8.2 includes a switch to plot ns ps calculated from n(x), p(x), jn(x),
% jp(x), jns and jps
% vX As with 8.2 but with more plot options
% TARR = array of time points
plotswitch = 1;     % 1 = plot n(x), p(x), 2 = plot ns, ps, 3 = 4x panels
calc_option = 2;    %1 = all terms of expression, 2 = first term only

par = sol.par;
T = par.T;
kB = par.kB;
dev = par.dev_ihalf;
u = sol.u;
t = sol.t;
x = sol.x;
x_sub = par.x_ihalf;
pcum1 = par.pcum + 1;   % Cumulative layer points array
al = par.active_layer;

xprime_p = dev.xprime_p;
xprime_n = dev.xprime_n;
xprime = dev.xprime;
alpha0 = dev.alpha0;
beta0 = dev.beta0;
mue = dev.mue;
muh = dev.muh;
dint = dev.dint;
Vt = par.kB*par.T;
kB = par.kB;
T = par.T;
pcum0 = par.pcum0;

dVdx = zeros(length(t), length(x_sub));
n = zeros(length(t), length(x_sub));
p = zeros(length(t), length(x_sub));
n_ana = zeros(length(t), length(x_sub));
p_ana = zeros(length(t), length(x_sub));
jn_ana = zeros(length(t), length(x_sub));
jp_ana = zeros(length(t), length(x_sub));

for i=1:length(t)
    [~, dVdx(i,:)] = pdeval(0, x, sol.u(i,:,1), x_sub);
    [n(i,:), ~] = pdeval(0, x, sol.u(i,:,2), x_sub);
    [p(i,:), ~] = pdeval(0, x, sol.u(i,:,3), x_sub);
end

alpha = par.q*dVdx./(par.kB*par.T) + alpha0;
beta = par.q*-dVdx./(par.kB*par.T) + beta0;

% Fluxes - approximate as on half mesh
[~, j, ~] = dfana.calcJ(sol);
jn = j.n;
jp = j.p;

% recombination
r_struct = dfana.calcr_ihalf(sol);
r = r_struct.tot;

%% Get indexes of interfaces
%int_index = find(contains(par.layer_type, 'interface'));   % only
%compatible from 2016 onwards
int_logical = strcmp(par.layer_type, 'interface');
loc = find(int_logical); % interface layer locations

ns = zeros(length(t), length(loc));     % Store time array of each ns in new column
ps = zeros(length(t), length(loc));     % Store time array of each ps in new column
jns = zeros(length(t), length(loc));     % Store time array of each ns in new column
jps = zeros(length(t), length(loc));     % Store time array of each ps in new column

for m = 1:length(loc)
    % Check location of interfacial surface carrier densities
    ns(:, m) = n(:, pcum1(loc(m)-1));
    jns(:, m) = jn(:, pcum1(loc(m)-1));     % jns and jps always from LH boundary
    ps(:, m) = p(:, pcum1(loc(m)-1));
    jps(:, m) = jp(:, pcum1(loc(m)-1));      % jns and jps always from LH boundary

    for kk = 1:length(tarr)
        k = find(sol.t <= tarr(kk));
        k = k(end);
        for i = 1:length(x_sub)
            if x_sub(i) >= par.dcum(loc(m)-1) && x_sub(i) <= par.dcum(loc(m))
                %% Carrier densities
                switch calc_option
                    case 1  % Full expressions
                        if alpha(k,i) <0
                            X = exp(alpha(k,i).*xprime(i));
                            n_ana(k,i) = ns(k,m)*X + jns(k,m)/(alpha(k,i)*mue(i)*Vt)*(1 - X)...
                                - ((jns(k,m)-jn(k,i))/(alpha(k,i)^2*mue(i)*Vt*xprime(i)))*(1 - X + alpha(k,i)*xprime(i));
                            
                        elseif alpha(k,i) > 0
                            Xstar = exp(alpha(k,i)*(xprime(i)-dint(k)));
                            n_ana(k,i) = ns(k,m)*Xstar + (jns(k,m)/(alpha(k,i)*mue(i)*Vt))*(1 - Xstar)...
                                - ((jns(k,m)-jn(k,i))/(alpha(k,i)^2*mue(i)*Vt*xprime(i)))*(1 + alpha(k,i)*xprime(i) - (1 + alpha(k,i)*dint(k))*Xstar);
                        end
                        
                        if beta(k,i) <0
                            Y = exp(beta(k,i)*xprime(i));
                            p_ana(k,i) = ps(k,m)*Y + (jps(k,m)/(beta(k,i)*muh(k)*Vt))*(1 - Y)...
                                - ((jps(k,m)-jp(k,i))/(beta(k,i)^2*muh(k)*Vt*xprime(i)))*(1 - Y + beta(k,i)*xprime(i));
                            
                        elseif beta(k,i) >= 0
                            Ystar = exp(beta(k,i)*(xprime(i)-dint(k)));
                            p_ana(k,i) = ps(k,m)*Ystar + (jps(k,m)/(beta(k,i)*muh(i)*Vt))*(1 - Ystar)...
                                - ((jps(k,m)-jp(k,i))/(beta(k,i)^2*muh(i)*Vt*xprime(i)))*(1 + beta(k,i)*xprime(i) - (1 + beta(k,i)*dint(k))*Ystar);
                        end
                        
                    case 2  % First terms only
                        if alpha(k,i) <0
                            X = exp(alpha(k,i).*xprime(i));
                            n_ana(k,i) = ns(k,m)*X;
                            
                        elseif alpha(k,i) > 0
                            Xstar = exp(alpha(k,i)*(xprime(i)-dint(k)));
                            n_ana(k,i) = ns(k,m)*Xstar;
                        end
                        
                        if beta(k,i) <0
                            Y = exp(beta(k,i)*xprime(i));
                            p_ana(k,i) = ps(k,m)*Y;
                            
                        elseif beta(k,i) >= 0
                            Ystar = exp(beta(k,i)*(xprime(i)-dint(k)));
                            p_ana(k,i) = ps(k,m)*Ystar;
                        end
                end
            end
            %% Fluxes
            jn_ana(k,i) = jns(k,m) - r(k,i)*(xprime(i));
            jp_ana(k,i) = jps(k,m) - r(k,i)*(xprime(i));        
        end 
    end
end

for m = 1:length(loc)
    if plotswitch == 1
        %% Plot carrier densities
        %% Interface 1
        figure(500)
        for i = 1:length(tarr)
            k = find(sol.t <= tarr(i));
            k = k(end);
            subplot(1,length(loc),m)
%             dfplot.colourblocks(sol, [1e-20,1e20]);
%             set(gca, 'YScale','log');
            semilogy(x_sub*1e7, n(k, :), x_sub*1e7, p(k, :),...
                x_sub*1e7, n_ana(k,:), 'k-.', x_sub*1e7, p_ana(k,:), 'k--')
            hold on
        end
        hold off
        title(['Interface ', num2str(m)])
        %legend('n', 'p', 'n-ana', 'p-ana')
        xlabel('Position, x (nm)')
        ylabel('Carrier density (cm-3)')
        xlim([1e7*par.dcum(loc(m)-1)-0.5,1e7*par.dcum(loc(m))+0.5])
        
        figure(600)
        %% Plot fluxes
        %% Interface 1
        for i = 1:length(tarr)
            k = find(sol.t <= tarr(i));
            k = k(end);
            subplot(1,length(loc),m)
            plot(x_sub*1e7, jn(k,:), x_sub*1e7, jp(k,:),...
                x_sub*1e7, jn_ana(i,:), 'k-.', x_sub*1e7, jp_ana(i,:), 'k--')
            hold on
        end
        title(['Interface ', num2str(m)])
        xlabel('Position, x (nm)')
        ylabel('Flux (cm-2s-1)')
        xlim([1e7*par.dcum(loc(m)-1)-0.5,1e7*par.dcum(loc(m))+0.5])
        hold off
        
    end
    
    if plotswitch ==2
        figure(700+m)
        for i = 1:length(tarr)
            k = find(sol.t <= tarr(i));
            k = k(end);
            subplot(1,length(loc),m)
            semilogy(x_sub*1e7, n(k, :), x_sub*1e7, p(k, :),...
                x_sub*1e7, n_ana(i,:), 'k-.', x_sub*1e7, p_ana(i,:), 'k--')
            hold on
        end
        hold off
        title(['Interface ', num2str(m)])
        %legend('n', 'p', 'n-ana', 'p-ana')
        xlabel('Position, x (nm)')
        ylabel('Carrier density (cm-3)')
        xlim([1e7*par.dcum(loc(m)-1)-0.5,1e7*par.dcum(loc(m))+0.5])
        
        for i = 1:length(tarr)
            k = find(sol.t <= tarr(i));
            k = k(end);
            subplot(1,length(loc),m)
            semilogy(x_sub*1e7, n(k, :), x_sub*1e7, p(k, :),...
                x_sub*1e7, n_ana(i,:), 'k-.', x_sub*1e7, p_ana(i,:), 'k--')
            hold on
        end
        hold off
        title(['Interface ', num2str(m)])
        %legend('n', 'p', 'n-ana', 'p-ana')
        xlabel('Position, x (nm)')
        ylabel('Carrier density (cm-3)')
        xlim([1e7*par.dcum(loc(m)-1)-0.5,1e7*par.dcum(loc(m))+0.5])
        
        % % components
        % figure(5011)
        % for i = 1:length(tarr)
        %     k = find(sol.t <= tarr(i));
        %     k = k(end);
        %     plot(x_sub*1e7, n1_comp1(i,:), x_sub*1e7, n1_comp2(i,:), x_sub*1e7, n1_comp3(i,:), x_sub*1e7, n1_ana(i,:), 'k-.')
        %     hold on
        % end
        % xlabel('Position, x (nm)')
        % ylabel('Carrier density (cm-3)')
        % legend('n comp', 'j comp', 'r comp', 'sum')
        % xlim([199.5,202.5])
        % hold off
    end
    
    if plotswitch == 3
        
        %% Plot carrier densities
        %% Interface 1
        figure(501)
        for i = 1:length(tarr)
            k = find(sol.t <= tarr(i));
            k = k(end);
            subplot(2,length(loc),m)
            semilogy(x_sub*1e7, n(k, :), x_sub*1e7, p(k, :),...
                x_sub*1e7, n_ana(i,:), 'k-.', x_sub*1e7, p_ana(i,:), 'k--')
            hold on
        end
        hold off
        title('Interface 1')
        %legend('n', 'p', 'n-ana', 'p-ana')
        xlabel('Position, x (nm)')
        ylabel('Carrier density (cm-3)')
        xlim([1e7*par.dcum(loc(m)-1)-0.5,1e7*par.dcum(loc(m))+0.5])
        
        %% Plot fluxes
        %% Interface 1
        for i = 1:length(tarr)
            k = find(sol.t <= tarr(i));
            k = k(end);
            subplot(2,length(loc),m)
            plot(x_sub*1e7, jn(k,:), x_sub*1e7, jp(k,:),...
                x_sub*1e7, jn_ana(i,:), 'k-.', x_sub*1e7, jp_ana(i,:), 'k--')
            hold on
        end
        title('Interface 1')
        xlabel('Position, x (nm)')
        ylabel('Flux (cm-2s-1)')
        xlim([1e7*par.dcum(loc(m)-1)-0.5,1e7*par.dcum(loc(m))+0.5])
        hold off
        
    end
end
end




