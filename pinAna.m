function [Voc, Vapp_arr, Jtot, Efn, Efp, U] = pinAna(solstruct)
% pinAna analyses the input solution and plots various useful graphs.
% Many plots are available to the user although currently these are
% commented out. In future

% Simple structure names
sol = solstruct.sol;
p = solstruct.p;

%% ANALYSIS %%
xnm = p.x*1e7; % x in nm for plotting

% split the solution into its component parts (e.g. electrons, holes and efield)
n = sol(:,:,1); % electrons
P = sol(:,:,2); % holes
a = sol(:,:,3); % mobile ions
V = sol(:,:,4); % electric potential

% Calculate energy levels and chemical potential
V = V - p.EA; % Electric potential
Ecb = p.EA-V-p.EA; % Conduction band potential
Evb = p.IP-V-p.EA; % Valence band potential
Efn = real(-V+p.Ei+(p.kB*p.T/p.q)*log(n/p.ni)); % Electron quasi-Fermi level
Efp = real(-V+p.Ei-(p.kB*p.T/p.q)*log(P/p.ni)); % Hole quasi-Fermi level

% interfacial points
p_i_array = ismembertol(p.x, p.tp);
i_n_array = ismembertol(p.x, p.tp + p.ti);

% indices of the interfacial points for delimiting integrals
p_i_index = find(p_i_array, 1);
i_n_index = find(i_n_array, 1);

if p.OC
    Voc = Efn(:, round(p.xpoints/2)) - Efp(:, 1); % Open Circuit Voltage
else
    Voc = p.Vapp;
end

%% TPV
if p.OC && p.pulseon % AC coupled mode
    Voc = Voc - Voc(1, :); % Removes baseline from TPV
    p.t = p.t-(p.pulsestart+p.pulselen); % Zero point adjustment
end

%% TPC
if ~p.OC && p.pulseon
    p.t = p.t-(p.pulsestart+p.pulselen); % TPC Zero point adjustment
end

%% Applied voltage array
switch p.JV % Current voltage array
    case 1
        Vapp_arr = p.Vstart + ((p.Vend-p.Vstart)*p.t*(1/p.tmax));
    case 2
        Vapp_arr = p.Vapp_func(p.Vapp_params, p.t);
    otherwise
        Vapp_arr = NaN;
end

%% Recombination
% separate charges by layer, the interfacial point will be present
% in more than one layer
n_h = n(:, 1:p_i_index);
P_h = P(:, 1:p_i_index);
n_i = n(:, p_i_index:i_n_index);
P_i = P(:, p_i_index:i_n_index);
n_e = n(:, i_n_index:end);
P_e = P(:, i_n_index:end);
U_h = trapz(p.x(1:p_i_index), p.kradhtl*(n_h.*P_h - p.ni^2) + (n_h.*P_h-p.ni^2)./(p.taun_htl.*(P_h+p.pthtl) + p.taup_htl.*(n_h+p.nthtl)), 2);
U_i = trapz(p.x(p_i_index:i_n_index), p.krad*(n_i.*P_i - p.ni^2) + (n_i.*P_i-p.ni^2)./(p.taun_i.*(P_i+p.pti) + p.taup_i.*(n_i+p.nti)), 2);
U_e = trapz(p.x(i_n_index:end), p.kradetl*(n_e.*P_e - p.ni^2) + (n_e.*P_e-p.ni^2)./(p.taun_etl.*(P_e+p.ptetl) + p.taup_etl.*(n_e+p.ntetl)), 2);
U = 1000*p.e*(U_h + U_i + U_e);

%% Current J

if ~p.OC
    % Calculates curkrent at every point and all times
    % Note the drift and diffusion currents do not cancel properly here
    switch p.calcJ
        case 0
            Jtot = NaN; % if calcJ is 0 the current is set to NaN
        case 1
            % find the internal current density in the device
            Jndiff = zeros(length(p.t), length(p.x));
            Jndrift = zeros(length(p.t), length(p.x));
            Jpdiff = zeros(length(p.t), length(p.x));
            Jpdrift = zeros(length(p.t), length(p.x));
            Jpart = zeros(length(p.t), length(p.x));
            Jidiff = zeros(length(p.t), length(p.x));
            Jidrift= zeros(length(p.t), length(p.x));
            dVdxt = zeros(length(p.t), length(p.x));

            for j=1:length(p.t)

                [nloc,dnlocdx] = pdeval(0,p.x,n(j,:),p.x);
                [Ploc,dPlocdx] = pdeval(0,p.x,P(j,:),p.x);
                [aloc,dalocdx] = pdeval(0,p.x,a(j,:),p.x);
                [~, dVdx] = pdeval(0,p.x,V(j,:),p.x);

                % Particle currents
                Jndiff(j,:) = (p.mue_i*p.kB*p.T*dnlocdx)*(1000*p.e);
                Jndrift(j,:) = (-p.mue_i*nloc.*dVdx)*(1000*p.e);

                Jpdiff(j,:) = (-p.muh_i*p.kB*p.T*dPlocdx)*(1000*p.e);
                Jpdrift(j,:) = (-p.muh_i*Ploc.*dVdx)*(1000*p.e);

                Jidiff(j,:) = (-p.mui*p.kB*p.T*dalocdx)*(1000*p.e);
                Jidrift(j,:) = (-p.mui*aloc.*dVdx)*(1000*p.e);

                % Particle current
                Jpart(j,:) = Jndiff(j,:) + Jndrift(j,:) + Jpdiff(j,:) + Jpdrift(j,:) + Jidiff(j,:) + Jidrift(j,:);

                % Electric Field
                dVdxt(j,:) = dVdx;
            end

            Jpartr = median(Jpart,2);

            % Displacement Current at right hand side
            Jdispr = (p.e*1000)*p.eppn*-gradient(dVdxt(:, end), p.t);

            % this is wrong
            %Jtot = Jpartr + Jdispr;
        case 2
            % Current calculation from the continuity equations
            [~, dndt] = gradient(n, p.x, p.t);
            dndt_tot = 1000*p.e*trapz(p.x, dndt, 2);

            %% Generation
            % Uniform Generation
            if p.OM == 0
                if p.Int ~= 0
                    g_tot = 1000*p.e*(p.x(i_n_index) - p.x(p_i_index)) * p.Int*p.G0;
                else
                    g_tot = 0;
                end
            end

            Jtot = dndt_tot - g_tot + U;

        case 3
            % like calcJ 1 but calculating just at right end boundary
        
            % find the internal current density in the device
            Jndiff = zeros(length(p.t), 1);
            Jndrift = zeros(length(p.t), 1);
            Jpdiff = zeros(length(p.t), 1);
            Jpdrift= zeros(length(p.t), 1);
            Jidiff = zeros(length(p.t), 1);
            Jidrift= zeros(length(p.t), 1);
            Jpart = zeros(length(p.t), 1);
            dVdxt = zeros(length(p.t), 1);

            for j=1:length(p.t)
                [nloc,dnlocdx] = pdeval(0,p.x,n(j,:),p.x(end));
                [Ploc,dPlocdx] = pdeval(0,p.x,P(j,:),p.x(end));
                [aloc,dalocdx] = pdeval(0,p.x,a(j,:),p.x(end));
                [~, dVdx] = pdeval(0,p.x,V(j,:),p.x(end));

                % Particle currents
                Jndiff(j) = (p.mue_n*p.kB*p.T*dnlocdx)*(1000*p.e);
                Jndrift(j) = (-p.mue_n*nloc.*dVdx)*(1000*p.e);

                Jpdiff(j) = (-p.muh_n*p.kB*p.T*dPlocdx)*(1000*p.e);
                Jpdrift(j) = (-p.muh_n*Ploc.*dVdx)*(1000*p.e);

                Jidiff(j) = (-p.mui*p.kB*p.T*dalocdx)*(1000*p.e);
                Jidrift(j) = (-p.mui*aloc.*dVdx)*(1000*p.e);

                % Particle current
                Jpart(j) = Jndiff(j) + Jndrift(j) + Jpdiff(j) + Jpdrift(j) + Jidiff(j) + Jidrift(j);

                % Electric Field
                dVdxt(j) = dVdx;
            end

            % Displacement Current at right hand side
            Jdispr = (p.e*1000)*p.eppn*gradient(dVdxt, p.t);

            Jtot = Jpart + Jdispr;
        
        case 4
            assert(p.BC == 3, [mfilename ' - calcJ 4 have to be used in combination with BC 3'])
            % Currents calculated from right-hand boundary
            Jn_r = p.sn_ext*(n(:, end) - p.etln0)*1000*-p.e;
            Jp_r = p.sp_rec*(P(:, end) - p.etlp0)*1000*p.e;
            Jn_l = p.sn_rec*(n(:, 1) - p.htln0)*1000*p.e;
            Jp_l = p.sp_ext*(P(:, 1) - p.htlp0)*1000*-p.e;

            Jtot = (Jn_r + Jp_r + Jn_l + Jp_l) / 2;

        case 5
            % like calcJ 1 but without the for cycle, faster. Gradient has a
            % slightly different output from pdeval.
            % Calculates curkrent at every point and all times
            % Note the drift and diffusion currents do not cancel properly here

            % p.t is not used but required
            dndx = gradient(n, p.x, p.t);
            dPdx = gradient(P, p.x, p.t);
            dadx = gradient(a, p.x, p.t);
            dVdx = gradient(V, p.x, p.t);

            % Particle currents
            Jndiff = (p.mue_i*p.kB*p.T*dndx)*(1000*p.e);
            Jndrift = (-p.mue_i*n.*dVdx)*(1000*p.e);

            Jpdiff = (-p.muh_i*p.kB*p.T*dPdx)*(1000*p.e);
            Jpdrift = (-p.muh_i*P.*dVdx)*(1000*p.e);

            Jidiff = (-p.mui*p.kB*p.T*dadx)*(1000*p.e);
            Jidrift = (-p.mui*a.*dVdx)*(1000*p.e);

            % Particle current
            Jpart = Jndiff + Jndrift + Jpdiff + Jpdrift + Jidiff + Jidrift;

            Jpartr = median(Jpart, 2);

            % Displacement Current at right hand side
            Jdispr = (p.e*1000)*p.eppn*-gradient(dVdx(:, end), p.t);

            Jtot = Jpartr + Jdispr;
        
        case 6
            % like calcJ 3 but without the for cycle, faster.

            % p.t is not used but required
            dndx_10p = gradient(n(:, end-10:end), p.x(end-10:end), p.t);
            dndx = dndx_10p(end);
            dPdx_10p = gradient(P(:, end-10:end), p.x(end-10:end), p.t);
            dPdx = dPdx_10p(end);
            dadx_10p = gradient(a(:, end-10:end), p.x(end-10:end), p.t);
            dadx = dadx_10p(end);
            dVdx_10p = gradient(V(:, end-10:end), p.x(end-10:end), p.t);
            dVdx = dVdx_10p(end);

            % Particle currents
            Jndiff = (p.mue_i*p.kB*p.T*dndx)*(1000*p.e);
            Jndrift = (-p.mue_i*n.*dVdx)*(1000*p.e);

            Jpdiff = (-p.muh_i*p.kB*p.T*dPdx)*(1000*p.e);
            Jpdrift = (-p.muh_i*P.*dVdx)*(1000*p.e);

            Jidiff = (-p.mui*p.kB*p.T*dadx)*(1000*p.e);
            Jidrift = (-p.mui*a.*dVdx)*(1000*p.e);

            % Particle current
            Jpart = Jndiff + Jndrift + Jpdiff + Jpdrift + Jidiff + Jidrift;

            % Displacement Current at right hand side
            Jdispr = (p.e*1000)*p.eppn*-gradient(dVdx, p.t);

            Jtot = Jpart + Jdispr;
        otherwise
            error([mfilename ' - calcJ value not recognized'])
    end
else
    Jtot = NaN; % if OC is 1 the current is set to NaN
end

%% GRAPHING %%

%% Spatial mesh
if p.meshx_figon == 1
    xmir = p.x;
    pxmir = 1:1:length(p.x);
    
    figure(1010);
        plot(xmir, pxmir, '.');
        xlabel('Position');
        ylabel('Point');
end

%% Time mesh
if p.mesht_figon == 1
    tmir = p.t;
    ptmir = 1:1:length(p.t);
    
    figure(200);
        plot(tmir, ptmir, '.');
        xlabel('Time');
        ylabel('Point');
end

%% General figures
if p.figson == 1
    
    % Open circuit voltage
    if p.OC == 1
        figure(7);
            plot (p.t, Voc);
            xlabel('Time [s]');
            ylabel('Voltage [V]');
    end
    
    % Defines end points for the graphing
    if p.OC
        xnmend = round(xnm(end)/2);
    else
        xnmend = xnm(end);
    end
    
    % Remove static ionic charge from contacts for plotting
    a_plotting = a(end, :) - p.NI * [ones(1, p_i_index-1), zeros(1, i_n_index - p_i_index + 1), ones(1, p.xpoints - i_n_index)];

    %% Band Diagram - subplot 1
    figure('Name', ['Band Diagram - ' inputname(1)], 'NumberTitle', 'off');
        subplot(3,1,1);
            plot (xnm, Efn(end,:), '--', xnm, Efp(end,:), '--', xnm, Ecb(end, :), xnm, Evb(end ,:));
            legend('E_{fn}', 'E_{fp}', 'CB', 'VB');
            set(legend,'FontSize',12);
            %xlabel('Position [nm]');
            ylabel('Energy [eV]');
            xlim([0, xnmend]);
            ylim([-3, 0.5]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            grid off;
    
        %% Electronic Charge Densities - subplot 2
        subplot(3,1,2); 
            semilogy(xnm, (sol(end,:,1)), xnm, (sol(end,:,2)));
            ylabel('{\itn, p} [cm^{-3}]')
            legend('\itn', '\itp')
            %xlabel('Position [nm]')
            xlim([0, xnmend]);
            ylim([1e0, 1e20]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            grid off
    
        %% Ionic charge density - subplot 3
        subplot(3,1,3);
            plot(xnm, a_plotting/1e19, 'black');
            ylabel('{\ita} [x10^{19} cm^{-3}]')
            xlabel('Position [nm]')
            xlim([0, xnmend]);
            ylim([0, 1.1*(max(sol(end,:,3))/1e19)]);
            %set(legend,'FontSize',12);
            %set(legend,'EdgeColor',[1 1 1]);
            grid off
    
    %% Currents as a function of position
    if p.calcJ == 1 || p.calcJ == 5
        figure(8);
            plot(xnm,Jndiff(end, :),xnm,Jndrift(end, :),xnm,Jpdiff(end, :),xnm,Jpdrift(end, :),xnm,Jidiff(end, :),xnm,Jidrift(end, :),xnm,Jpart(end, :));
            legend('Jn diff','Jn drift','Jp diff','Jp drift','Ji diff','Ji drift','Total J');
            xlabel('Position [nm]');
            ylabel('Current Density [mA cm^-2]');
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            xlim([0, xnmend]);
            grid off;
            drawnow;
    end
    
    %% Currents as a function of time
    if p.calcJ
        % Particle currents as a function of time
        figure(10);
            plot(p.t, Jtot);
            legend('J_n')
            xlabel('time [s]');
            ylabel('J [mA cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
            grid off;
            drawnow;
        
        if p.JV == 1
            figure(11)
                plot(Vapp_arr, Jtot)
                xlabel('V_{app} [V]')
                ylabel('Current Density [mA cm^-2]');
                grid off;
        end
    end
    drawnow
end

end

