function [Voc, Vapp_arr, Jn, Efn, Efp] = pinAna(solstruct)

% pinAna analyses the input solution and plots various useful graphs.
% Many plots are available to the user although currently these are
% commented out. In future 

% Simple structure names
sol = solstruct.sol;
p = solstruct.p;

%% ANALYSIS %%
xnm = p.x*1e7;    % x in nm for plotting

% split the solution into its component parts (e.g. electrons, holes and efield)
n = sol(:,:,1);     % electrons
P = sol(:,:,2);     % holes
a = sol(:,:,3);     % mobile ions
V = sol(:,:,4);     % electric potential

% Calculate energy levels and chemical potential         
V = V - p.EA;                                % Electric potential
Ecb = p.EA-V-p.EA;                             % Conduction band potential
Evb = p.IP-V-p.EA;                             % Valence band potential
Efn = real(-V+p.Ei+(p.kB*p.T/p.q)*log(n/p.ni));      % Electron quasi-Fermi level 
Efp = real(-V+p.Ei-(p.kB*p.T/p.q)*log(P/p.ni));      % Hole quasi-Fermi level
Phin = real(p.Ei+(p.kB*p.T/p.q)*log(n/p.ni)-p.EA);     % Chemical Potential electron
Phip = real(p.Ei-(p.kB*p.T/p.q)*log(P/p.ni)-p.EA);
Phi = Phin - Phip;

% p-type binary matrix
pBM = ones(length(p.t), p.xpoints)*diag(p.x < p.tp);
% Intrinsic binary matrix
iBM = ones(length(p.t), p.xpoints)*diag(p.x >= p.tp & p.x <= (p.tp + p.ti));
% n-type binary matrix
nBM = ones(length(p.t), p.xpoints)*diag(p.x > p.tp + p.ti);

nstat = zeros(1, p.xpoints);                                  % Static charge array
nstat = (-p.NA-p.NI)*pBM + (-p.NI*iBM) + (p.ND-p.NI)*nBM; %(-NA+pthtl-nthtl-NI)*pBM + (-NI*iBM) + (ND+ptetl-ntetl-NI)*nBM;   
rhoc = (-n + P + a + nstat);     % Net charge density calculated from adding individual charge densities

% Remove static ionic charge from contacts for plotting
a = a - (p.NI*pBM + p.NI*nBM);

% Recomination Rate - NEEDS SORTING 24/03/2016
Urec = 0;% (krad((n*p)-ni^2) + sn*(n-htln0)).*nBM + (krad((n*p)-ni^2)+ sp*(p-etlp0)).*pBM;

if p.OC == 1
    
    Voc = Efn(:, round(p.xpoints/2)) - Efp(:, 1);                    % Open Circuit Voltage
    Voc_chem = Phin(:, round(p.xpoints/2)) - Phip(:, 1);              % Chemical componenet
    Voc_V = V(:, round(p.xpoints/2)) - V(:, 1);

else
    
    Voc = nan;
    
end

%% TPV
if p.OC == 1  && p.pulseon == 1                 % AC coupled mode
   
    Voc = Voc - Voc(1, :);                  % Removes baseline from TPV
    p.t = p.t-(p.pulsestart+p.pulselen);            % Zero point adjustment                               
end

%% TPC
if p.OC == 0 && p.pulseon == 1 

    p.t = p.t-(p.pulsestart+p.pulselen);                    % TPC Zero point adjustment   

end

Potp = V(end, :);                                   % Potential 

rhoctot = trapz(p.x, rhoc, 2)/p.xmax;       % Integrated net charge

rho_a = a - p.NI;                         % Net ionic charge
rho_a_tot = trapz(p.x, rho_a, 2)/p.xmax;    % Total Net ion charge

ntot = trapz(p.x, n, 2);                  % Integrated electron density 
ptot = trapz(p.x, P, 2);                  % Integrated hole density

if p.JV == 1
    Vapp_arr = p.Vstart + ((p.Vend-p.Vstart)*p.t*(1/p.tmax)); % Current voltage array
else
    Vapp_arr = nan;
end

if p.OC == 0
    %% Current calculation from the continuity equations
    for j = 1:size(n, 2)
        
        dndt(:,j) = gradient(n(:,j), p.t);
        dpdt(:,j) = gradient(P(:,j), p.t);
        
    end
    
    dndtInt = trapz(p.x, dndt, 2);
    dpdtInt = trapz(p.x, dpdt, 2);
    
    %% Recombination
    Ubtb = p.krad*(n.*P - p.ni^2);
    Usrh = ((n.*P-p.ni^2)./((p.taun_htl.*(P+p.pthtl)) + (p.taup_htl.*(n+p.nthtl)))).*pBM...
        +((n.*P-p.ni^2)./((p.taun_etl.*(P+p.ptetl)) + (p.taup_etl.*(n+p.ntetl)))).*nBM;
    
    U = Ubtb + Usrh;
    
    %% Generation
    % Uniform Generation
    if p.OM == 0
        
        if p.Int ~= 0
            
            g = p.Int*p.G0*iBM;
            
        else
            
            g = 0;
            
        end
        
    end
    
    dJndx = dndt - g + U;
    
    Jn = trapz(p.x, dJndx, 2)*1000*p.e;
    
    %% Calculates current at every point and all times
    % Note the drift and diffusion currents do not cancel properly here
    if p.calcJ == 1
        
        % find the internal current density in the device
        Jndiff = zeros(length(p.t), length(p.x));
        Jndrift = zeros(length(p.t), length(p.x));
        Jpdiff = zeros(length(p.t), length(p.x));
        Jpdrift = zeros(length(p.t), length(p.x));
        Jpart = zeros(length(p.t), length(p.x));
        Jtot = zeros(length(p.t));
        
        for j=1:length(p.t)
            
            %tj = t(j);
            
            [nloc,dnlocdx] = pdeval(0,p.x,n(j,:),p.x);
            [ploc,dplocdx] = pdeval(0,p.x,P(j,:),p.x);
            [iloc,dilocdx] = pdeval(0,p.x,a(j,:),p.x);
            [Vloc, dVdx] = pdeval(0,p.x,V(j,:),p.x);
            
            % Particle currents
            Jndiff(j,:) = (p.mue_i*p.kB*p.T*dnlocdx)*(1000*p.q);
            Jndrift(j,:) = (-p.mue_i*nloc.*dVdx)*(1000*p.q);
            
            Jpdiff(j,:) = (-p.muh_i*p.kB*p.T*dplocdx)*(1000*p.q);
            Jpdrift(j,:) = (-p.muh_i*ploc.*dVdx)*(1000*p.q);
            
            Jidiff(j,:) = (-p.mui*p.kB*p.T*dilocdx)*(1000*p.q);
            Jidrift(j,:) = (-p.mui*iloc.*dVdx)*(1000*p.q);
            
            % Particle current
            Jpart(j,:) = Jndiff(j,:) + Jndrift(j,:) + Jpdiff(j,:) + Jpdrift(j,:) + Jidiff(j,:) + Jidrift(j,:);
            
            % Electric Field
            dVdxt(j,:) = dVdx;
            
        end
        
    end

else
    
    Jn = 0;
    
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
if p.OC == 1
    
    xnmend = round(xnm(end)/2);
    
else
    
    xnmend = xnm(end);
end

%% Band Diagram - subplot 1
FH1 = figure(1);
%set(FigHandle, 'units','normalized','position',[.1 .1 .4 .4]);
PH1 = subplot(3,1,1);
plot (xnm, Efn(end,:), '--', xnm, Efp(end,:), '--', xnm, Ecb(end, :), xnm, Evb(end ,:));
%legend('E_{fn}', 'E_{fp}', 'CB', 'VB');
set(legend,'FontSize',12);
%xlabel('Position [nm]');
ylabel('Energy [eV]'); 
xlim([0, xnmend]);
ylim([-3, 0.5]);
set(legend,'FontSize',12);
set(legend,'EdgeColor',[1 1 1]);
grid off;

%% Electronic Charge Densities - subplot 2
PH2 = subplot(3,1,2);
semilogy(xnm, (sol(end,:,1)), xnm, (sol(end,:,2)));
ylabel('{\itn, p} [cm^{-3}]')
%legend('\itn', '\itp')
%xlabel('Position [nm]')
xlim([0, xnmend]);
ylim([1e0, 1e20]);
set(legend,'FontSize',12);
set(legend,'EdgeColor',[1 1 1]);
grid off

%% Ionic charge density - subplot 3
PH3 = subplot(3,1,3);
plot(xnm, a(end,:)/1e19, 'black');
ylabel('{\ita} [x10^{19} cm^{-3}]')
xlabel('Position [nm]')
xlim([0, xnmend]);
ylim([0, 1.1*(max(sol(end,:,3))/1e19)]);
set(legend,'FontSize',12);
set(legend,'EdgeColor',[1 1 1]);
grid off

%% Ionic charge density
%{
figure(3)
%set(FigHandle, 'units','normalized','position',[.1 .1 .4 .4]);
plot(xnm, Irho (end,:))
legend('Ions')
set(legend,'FontSize',16);
xlabel('Position [nm]');
ylabel('Density [cm^{-3}]');
xlim([0, xnmend]);
set(legend,'FontSize',14);
set(legend,'EdgeColor',[1 1 1]);
grid off;
drawnow;
%}

%% Electric field potential surface
%{
figure(4);
surf(x,t,V);
title('Electric Potential (x,t)');
xlabel('Distance x');
ylabel('time');
%}

%% Space charge
%{
figure(5)
plot(xnm, rhoc(end, :))
ylabel('Space Charge Density [cm^{-3}]')
xlabel('Position [nm]')
xlim([0, xnmend]);
set(legend,'FontSize',14);
set(legend,'EdgeColor',[1 1 1]);
grid off
%}

%% Electron and hole charge densities
%{
figure(6)
plot(t, ntot, t, ptot)
ylabel('Charge Density [cm^{-3}]')
xlabel('time [s]')
legend('electrons', 'holes')
set(legend,'FontSize',14);
set(legend,'EdgeColor',[1 1 1]);
grid off
%}

%% Generation profile- not currently implemented
if p.OM == 1 && p.Int~=0 || p.OM == 2 && p.Int~=0

    genspacenm = genspace * 1e7;

    figure(7);
    plot(genspacenm, Gx1S, genspacenm, GxLas)
    ylabel('Generation Rate [cm^{3}s^{-1}]');
    xlabel('Position [nm]');
    legend('1 Sun', '638 nm');
    xlim([0, genspacenm(end)]);
    grid off

end

%% Currents as a function of position
if p.calcJ == 1

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

%{
% Electric Field
figure(9);
surf(xnm, t, dVdxt);
xlabel('Position [m]');
ylabel('time [s]');
title('Electric Field');
%}

end

%% Currents as a function of time
if p.calcJ == 0 || p.calcJ == 1

% Particle and displacement currents as a function of time
figure(10);
plot(p.t, Jn);
legend('J_n')%, 'Jparticle', 'Jdisp')
xlabel('time [s]');
ylabel('J [mA cm^{-2}]');
set(legend,'FontSize',16);
set(legend,'EdgeColor',[1 1 1]);
grid off;
drawnow;

    if p.JV == 1

        figure(11)
        plot(Vapp_arr, Jn)
        xlabel('V_{app} [V]')
        ylabel('Current Density [mA cm^-2]');
        grid off;
    
    end

end

drawnow

end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unused figures

%{
figure(1);
surf(x,t,n);
set(gca, 'ZScale', 'log')
%semilogy(x,n(end,:));
title('n(x,t)');
xlabel('Distance x');
ylabel('time');

figure(2);
surf(x,t,p);
set(gca, 'ZScale', 'log')
title('p(x,t)');
xlabel('Distance x');
ylabel('time');

figure(3);
surf(x,t,a);
title('ion(x,t)');
xlabel('Distance x');
ylabel('time');

figure(11)
plot(xnm, Fc, xnm, Fp);
xlim([0, (xnmend/2)]);
legend('from charge', 'from pot')
ylabel('E Field [V/cm]');
grid off;

% Electric Field vs Position
figure(6);
plot(xnm, Fp(end, :));
xlabel('Position [nm]');
ylabel('Electric Field [Vcm^{-1}]');
grid off;

%}

%{
figure(3)
%set(FigHandle, 'units','normalized','position',[.1 .1 .4 .4]);
[Ax1, h1, h2] = plotyy(xnm, rhoc(end, :), xnm, Urec(end, :));
linkaxes(Ax1,'x');
set(gca,'xlim',[0, (xnmend/2)]);
ylabel(Ax1(1),'Net Charge Density [cm^{-3}]') % left y-axis
ylabel(Ax1(2),{'Recombination';'Rate [cm^{-3}s^{-1}]'}) % right y-axis
%set(Ax1(2),'YScale','log')
set(Ax1(1),'Position', [0.1 0.11 0.7 0.8]);
set(Ax1(2),'Position', [0.1 0.11 0.7 0.8]);
set(Ax1(1),'ycolor',[0.1, 0.1, 0.1]) 
set(Ax1(2),'ycolor',[0, 0.4470, 0.7410])
set(h1,'color',[0.1, 0.1, 0.1])
set(h2, 'color',[0, 0.4470, 0.7410])
grid off;

figure(4)
plot(xnm, g);
grid off

% Electric and Checmial Potential components
figure(4)
%set(FigHandle, 'units','normalized','position',[.1 .1 .4 .4]);
[Ax2, h3, h4] = plotyy(xnm, Potp, xnm, Phi(end,:));
linkaxes(Ax2,'x');
set(gca,'xlim',[0, (xnmend/2)]);
ylabel(Ax2(1),'Electric Potential [V]') % left y-axis
ylabel(Ax2(2),'Chemical Potential [V]') % right y-axis
% get current (active) axes property
set(Ax2(1),'Position', [0.13 0.11 0.775-.08 0.815]);
set(Ax2(2),'Position', [0.13 0.11 0.775-.08 0.815]);
set(Ax2(1),'ycolor',[0.1, 0.1, 0.1]) 
set(Ax2(2),'ycolor',[0.8500, 0.3250, 0.0980])
set(h3,'color',[0.1, 0.1, 0.1])
grid off;



% figure(200)
% [AX,H1,H2] = plotyy(xnm, [Jidiff(end, :).', Jidrift(end, :).'], xnm, (sol(end,:,3)));
% legend('Ion diffusion', 'Ion drift', 'a')
% xlabel('Position [nm]');
% ylabel('Density [cm^{-3}]/Current Density');
% set(AX(2),'Yscale','linear');
% set(legend,'FontSize',20);
% set(legend,'EdgeColor',[1 1 1]);
% set(AX(1), 'Position',[0.18 0.18 0.7 0.70]);     % left, bottom, width, height       
% set(AX(2), 'Position',[0.18 0.18 0.7 0.70]);
% box on
% set(AX(1), 'YMinorTick','on');     % left, bottom, width, height       
% set(AX(2), 'XMinorTick','on','YMinorTick','on');
% set(AX(1),'xlim',[190 250]);
% set(AX(2),'xlim',[190 250]);
% %set(AX(1),'ylim',[1e6 1e18]);
% set(AX(2),'ycolor',[0.9290    0.6940    0.1250]);
% set(H2,'color',[0.9290    0.6940    0.1250]);
% set(legend,'FontSize',12);
% set(legend,'EdgeColor',[1 1 1]);
% grid off;

%}

%% Delete unwanted fields from params
