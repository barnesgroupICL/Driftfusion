function [Voc, Vapp_arr, Jtotr] = pinAna(solstruct)

% pinAna analyses the input solution and plots various useful graphs.
% Many plots are available to the user although currently these are
% commented out. In future 

% Simple structure names
sol = solstruct.sol;
params = solstruct.params;

% Unpacks params structure for use in current workspace 
v2struct(params);

%% ANALYSIS %%
xnm = x*1e7;    % x in nm for plotting

% split the solution into its component parts (e.g. electrons, holes and efield)
n = sol(:,:,1);     % electrons
p = sol(:,:,2);     % holes
a = sol(:,:,3);     % mobile ions
V = sol(:,:,4);     % electric potential

% Calculate energy levels and chemical potential         
V = V - EA;                                % Electric potential
Ecb = EA-V-EA;                             % Conduction band potential
Evb = IP-V-EA;                             % Valence band potential
Efn = real(-V+Ei+(kB*T/q)*log(n/ni));      % Electron quasi-Fermi level 
Efp = real(-V+Ei-(kB*T/q)*log(p/ni));      % Hole quasi-Fermi level
Phin = real(Ei+(kB*T/q)*log(n/ni)-EA);     % Chemical Potential electron
Phip = real(Ei-(kB*T/q)*log(p/ni)-EA);
Phi = Phin - Phip;

% p-type binary matrix
pBM = ones(length(t), xpoints)*diag(x < tp);
% Intrinsic binary matrix
iBM = ones(length(t), xpoints)*diag(x >= tp & x <= (tp + ti));
% n-type binary matrix
nBM = ones(length(t), xpoints)*diag(x > tp + ti);

nstat = zeros(1, xpoints);                                  % Static charge array
nstat = (-NA-NI)*pBM + (-NI*iBM) + (ND-NI)*nBM; %(-NA+pthtl-nthtl-NI)*pBM + (-NI*iBM) + (ND+ptetl-ntetl-NI)*nBM;   
rhoc = (-n + p + a + nstat);     % Net charge density calculated from adding individual charge densities

% Remove static ionic charge from contacts for plotting
a = a - (NI*pBM + NI*nBM);

% Recomination Rate - NEEDS SORTING 24/03/2016
Urec = 0;% (krad((n*p)-ni^2) + sn*(n-htln0)).*nBM + (krad((n*p)-ni^2)+ sp*(p-etlp0)).*pBM;

if OC == 1
    
    Voc = Efn(:, round(xpoints/2)) - Efp(:, 1);                    % Open Circuit Voltage
    Voc_chem = Phin(:, round(xpoints/2)) - Phip(:, 1);              % Chemical componenet
    Voc_V = V(:, round(xpoints/2)) - V(:, 1);

else
    
    Voc = nan;
    
end

%% TPV
if OC == 1  && pulseon == 1                 % AC coupled mode
   
    Voc = Voc - Voc(1, :);                  % Removes baseline from TPV
    t = t-(pulsestart+pulselen);            % Zero point adjustment                               
end

%% TPC
if OC == 0 && pulseon == 1 

    t = t-(pulsestart+pulselen);                    % TPC Zero point adjustment   

end

Potp = V(end, :);                                   % Potential 

rhoctot = trapz(x, rhoc, 2)/xmax;       % Integrated net charge

rho_a = a - NI;                         % Net ionic charge
rho_a_tot = trapz(x, rho_a, 2)/xmax;    % Total Net ion charge

ntot = trapz(x, n, 2);                  % Integrated electron density 
ptot = trapz(x, p, 2);                  % Integrated hole density

if JV == 0
    
    Vapp_arr = nan;                     % Current voltage array
    
elseif JV == 1
    
    Vapp_arr = Vstart + ((Vend-Vstart)*t*(1/tmax));
    
end

if calcJ == 0
    
    Jtotr = nan;
    

% Calculates current at every point and all times
elseif calcJ == 1

% find the internal current density in the device
Jndiff = zeros(length(t), length(x));
Jndrift = zeros(length(t), length(x));
Jpdiff = zeros(length(t), length(x));
Jpdrift = zeros(length(t), length(x));
Jpart = zeros(length(t), length(x));
Jtot = zeros(length(t));   
    
for j=1:length(t)
    
    %tj = t(j);
    
    [nloc,dnlocdx] = pdeval(0,x,n(j,:),x);    
    [ploc,dplocdx] = pdeval(0,x,p(j,:),x);
    [iloc,dilocdx] = pdeval(0,x,a(j,:),x);
    [Vloc, dVdx] = pdeval(0,x,V(j,:),x);
    
    % Particle currents
    Jndiff(j,:) = (mue_i*kB*T*dnlocdx)*(1000*e);
    Jndrift(j,:) = (-mue_i*nloc.*dVdx)*(1000*e);
   
    Jpdiff(j,:) = (-muh_i*kB*T*dplocdx)*(1000*e);
    Jpdrift(j,:) = (-muh_i*ploc.*dVdx)*(1000*e);
    
    Jidiff(j,:) = (-mui*kB*T*dilocdx)*(1000*e);
    Jidrift(j,:) = (-mui*iloc.*dVdx)*(1000*e);

    % Particle current
    Jpart(j,:) = Jndiff(j,:) + Jndrift(j,:) + Jpdiff(j,:) + Jpdrift(j,:) + Jidiff(j,:) + Jidrift(j,:);   
    
    % Electric Field
    dVdxt(j,:) = dVdx;
    
end

% Currents at the boundaries (should be the same)
%Jpartl = Jpart(:,1);
%Jpartr = Jpart(:,end);
Jpartr = median(Jpart,2);
Jpartr = Jpartr'; 

% Displacement Current at right hand side
Fend = -(dVdxt(:, end));
Jdispr = (e*1000)*eppn*-gradient(dVdxt(:, end), t);
Jdispr = Jdispr';

Jtotr = Jpartr + Jdispr;    

% Calculates currents only for right hand x points at all times
elseif calcJ == 2
    
% find the internal current density in the device
Jndiff = zeros(length(t), 1);
Jndrift = zeros(length(t), 1);
Jpdiff = zeros(length(t), 1);
Jpdrift= zeros(length(t), 1);
Jpart = zeros(length(t), 1);
    
    for j=1:length(t)
        
    [nloc,dnlocdx] = pdeval(0,x,n(j,:),x(end));    
    [ploc,dplocdx] = pdeval(0,x,p(j,:),x(end));
    [iloc,dilocdx] = pdeval(0,x,a(j,:),x(end));
    [Vloc, dVdx] = pdeval(0,x,V(j,:),x(end));
    
    % Particle currents
    Jndiff(j) = (mue_n*kB*T*dnlocdx)*(1000*e);
    Jndrift(j) = (-mue_n*nloc.*dVdx)*(1000*e);
   
    Jpdiff(j) = (-muh_n*kB*T*dplocdx)*(1000*e);
    Jpdrift(j) = (-muh_n*ploc.*dVdx)*(1000*e);
    
    Jidiff(j) = (-mui*kB*T*dilocdx)*(1000*e);
    Jidrift(j) = (-mui*iloc.*dVdx)*(1000*e);

    % Particle current
    Jpart(j) = Jndiff(j) + Jndrift(j) + Jpdiff(j) + Jpdrift(j) + Jidiff(j) + Jidrift(j);   
    
    % Electric Field
    dVdxt(j) = dVdx;

    end

% Currents at the boundaries
Jpartr = Jpart';

%Jpartr = -(sn*n(:, end) - ni) %Check when surface recombination is used

% Displacement Current at right hand side

Jdispr = (e*1000)*eppn*gradient(dVdxt, t);

Jtotr = Jpartr + Jdispr;     

if pulseon == 1
    
    Jtotr = Jtotr - Jtotr(end);    % remove baseline

end

% Current calculated from QFL
elseif calcJ == 3

        for j=1:length(t)

            dEfndx(j,:) = gradient(Efn(j, :), x);
            dEfpdx(j,:) = gradient(Efp(j, :), x);

            [Vloc, dVdx] = pdeval(0,x,V(j,:),x(end));
             % Electric Field
            dVdxt(j) = dVdx;

        end

    Jpart =  mue_i*n.*dEfndx*(1000*e) +  muh_i*p.*dEfpdx*(1000*e);

    Jdispr = (e*1000)*eppn*gradient(dVdxt, t);
    Jpartr = Jpart(:,pe+0.2*pp);
    Jtotr = Jpartr + Jdispr;
    Jdispr = 0;

end

%% GRAPHING %%

%% Spatial mesh
if meshx_figon == 1
    
    xmir = x;
    pxmir = 1:1:length(x);
    
    figure(1010);
    plot(xmir, pxmir, '.');
    xlabel('Position');
    ylabel('Point');
    
end


%% Time mesh
if mesht_figon == 1

    tmir = t;
    ptmir = 1:1:length(t);
    
    figure(200);
    plot(tmir, ptmir, '.');
    xlabel('Time');
    ylabel('Point');

end

%% General figures
if figson == 1;
    
    % Open circuit voltage
      if OC == 1
        
        figure(7);
        plot (t, Voc);
        xlabel('Time [s]');   
        ylabel('Voltage [V]');

      end

% Defines end points for the graphing
if OC == 1
    
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
if OM == 1 && Int~=0 || OM == 2 && Int~=0

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
if calcJ == 1

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
if calcJ == 1 || calcJ == 2 || calcJ == 3

% Particle and displacement currents as a function of time
figure(10);
plot(t, Jtotr, t, Jpartr, t, Jdispr);
legend('Jtotal', 'Jparticle', 'Jdisp')
xlabel('time [s]');
ylabel('J [mA cm^{-2}]');
set(legend,'FontSize',16);
set(legend,'EdgeColor',[1 1 1]);
grid off;
drawnow;

    if JV == 1
        %JV
        figure(11)
        plot(Vapp_arr, Jtotr)
        xlabel('V_{app} [V]')
        ylabel('Current Density [mA cm^-2]');
        grid off;
    end

end

drawnow

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




end