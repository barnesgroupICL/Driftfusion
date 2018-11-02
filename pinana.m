function [Voc, Vapp_arr, Jtot] = pinana(varargin)

% Plotting defaults
set(0,'DefaultLineLinewidth',1);
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
set(0,'DefaultAxesXcolor', [0, 0, 0]);
set(0,'DefaultAxesYcolor', [0, 0, 0]);
set(0,'DefaultAxesZcolor', [0, 0, 0]);
set(0,'DefaultTextColor', [0, 0, 0]);

% tarr is a time time array for the time you wish to plot
if length(varargin) == 1
    solstruct = varargin{1};
    tarr = solstruct.t(end);
    pointtype = 't';
elseif length(varargin) == 2
    solstruct = varargin{1};
    tarr = varargin{2};
    pointtype = 't';
elseif length(varargin) == 3
    solstruct = varargin{1};
    pointtype = varargin{2};
    tarr = varargin{3};
end

% Simple structure names
sol = solstruct.sol;
par = solstruct.p;
x = solstruct.x;
t = solstruct.t;

if par.OM == 1
    gx = solstruct.gx;
end

par.x = x;
par.t = t;        % For backwards compatibility

xpoints = length(x);

%% ANALYSIS %%
xnm = x*1e7;    % x in nm for plotting

%%%%% ANALYSIS %%%%%

% split the solution into its component parts (e.g. electrons, holes and efield)
n = sol(:,:,1);
p = sol(:,:,2);
a = sol(:,:,3);
V = sol(:,:,4);

% Create 2D matrices for multiplication with solutions
EAmat = repmat(par.dev.EA, length(t), 1);
IPmat = repmat(par.dev.IP, length(t), 1);
muemat = repmat(par.dev.mue, length(t), 1);
muhmat = repmat(par.dev.muh, length(t), 1);
muionmat = repmat(par.dev.muion, length(t), 1);
NAmat = repmat(par.dev.NA, length(t), 1);
NDmat = repmat(par.dev.ND, length(t), 1);
N0mat = repmat(par.dev.N0, length(t), 1);
Nionmat = repmat(par.dev.Nion, length(t), 1);
eppmat = repmat(par.dev.epp, length(t), 1);
nimat = repmat(par.dev.ni, length(t), 1);
kradmat = repmat(par.dev.krad, length(t), 1);

Ecb = EAmat-V;                                 % Conduction band potential
Evb = IPmat-V;                                 % Valence band potential

Efn = zeros(size(n,1), size(n,2));
Efp = zeros(size(n,1), size(n,2));

if par.stats == 'Fermi'
    
    for i = 1:size(n,1)           % time
        for j = 1:size(n,2)       % position
            Efn(i,j) = F.Efn_fd_fun(n(i,j), par.dev.Efn(j,:),  par.dev.n_fd(j,:));
            Efp(i,j) = F.Efp_fd_fun(p(i,j), par.dev.Efp(j,:),  par.dev.p_fd(j,:));
        end
    end
    Efn = Efn-V;
    Efp = Efp-V;
    
elseif par.stats == 'Boltz'
    Efn = real(Ecb+(par.kB*par.T/par.q)*log(n./N0mat));        % Electron quasi-Fermi level
    Efp = real(Evb-(par.kB*par.T/par.q)*log(p./N0mat));        % Hole quasi-Fermi level
end

% Remove ionic charge densities from contact regions
rhoa = a - Nionmat;

%nstat = zeros(1, xpoints);                                  % Static charge array
rhostat = NAmat+NDmat;
rhoc = (-n + p + rhostat);     % Net charge density calculated from adding individual charge densities

Voc = nan;

if par.OC == 1  && par.pulseon == 1                               % AC coupled mode
    
    Voc = Voc - Voc(1, :);
    t = (t-(par.pulsestart+par.pulselen));          % Zero point adjustment
end

if par.OC == 0 && par.pulseon == 1
    
    t = (t-par.pulsestart);          % Zero point adjustment
    
end


for i=1:length(t)
    
    Fp(i,:) = -gradient(V(i, :), x);                      % Electric field calculated from V
    
end

Potp = V(end, :);

rhoctot = trapz(x, rhoc, 2)/par.dcum(end);   % Net charge

Irho = a - Nionmat;                  % Net ionic charge
Irhotot = trapz(x, Irho, 2)/par.dcum(end);   % Total Net ion charge

ntot = trapz(x, n, 2);     % Total
ptot = trapz(x, p, 2);

if par.JV == 1
    
    Vapp_arr = par.Vstart + ((par.Vend-par.Vstart)*t*(1/par.tmax));
    
else
    
    Vapp_arr = nan;
    
end

%% Current calculation from continuity equations

for j = 1:size(n, 2)
    
    dndt(:,j) = gradient(n(:,j), t);
    dpdt(:,j) = gradient(p(:,j), t);
    
end

dndtInt = trapz(x, dndt, 2);
dpdtInt = trapz(x, dpdt, 2);

% Recombination
Ubtb = kradmat.*(n.*p - nimat.^2);

Usrh = 0;%((n.*p - par.ni(2)^2)./((par.taun(1).*(p+par.pt(2))) + (par.taup(1).*(n+par.nt(2))))).*piBM...
% + ((n.*p- par.ni(2)^2)./((par.taun(3).*(p+par.pt(2))) + (par.taup(3).*(n+par.nt(2))))).*inBM;

U = Ubtb + Usrh;

% Uniform Generation
switch par.OM
    
    % Uniform generation
    case 0
        
        g = par.Int*par.dev.G0;
        
    case 1
        
        gxAM15 = par.Int*repmat(gx.AM15', length(t), 1);
        
        if par.pulseon == 1
            
            las = repmat(gx.las', length(t), 1);
            pulset = ones(length(x), length(t))*diag(t >= par.pulsestart & t < par.pulselen + par.pulsestart);
            pulset = pulset';
            gxlas = las.*pulset;
            
        else
            gxlas = 0;
            
        end
        
        g = gxAM15 + gxlas;
        
    case 2
        % Transfer Matrix
        if par.Int == 0
            
            g = 0;
            
        else
            
            g = par.Int*interp1(par.genspace, solstruct.Gx1S, (x-par.dcum(1)));
            
        end
        
end

djndx = -(dndt - g + U);    % Not certain about the sign here
djpdx = -(dpdt - g + U);

% Integrate across the device to get delta fluxes at all positions
deltajn = cumtrapz(par.x, djndx, 2);
deltajp = cumtrapz(par.x, djpdx, 2);

%% Currents from the boundaries
    switch par.BC
        case 0
            jn_l = 0;
            jp_l = 0;
            jn_r = 0;
            jp_r = 0;
            % Blocking contacts
        case 1
            % Setting jp_l = djpdx(end) ensures that jp_r = 0;
            jn_l = 0;
            jp_l = -deltajp(:, end);
            
            jn_r = deltajn(:, end);
            jp_r = 0;
            
        case 2
            
            jn_l = -par.sn_l*(n(:, 1) - par.nleft);
            jp_l = -deltajp(:, end) + par.sp_r*(p(:, end) - par.pright);
            
            jn_r = deltajn(:, end) - par.sn_l*(n(:, 1) - par.nleft);
            jp_r = par.sp_r*(p(:, end) - par.pright);
            
        case 3
            
            jn_l = -par.sn_l*(n(:, 1) - par.nleft);
            jp_l = -par.sp_l*(p(:, 1) - par.pleft);
            
            jn_r = par.sn_r*(n(:, end) - par.nright);
            jp_r = par.sp_r*(p(:, end) - par.pright);
            
    end

% Calculate total electron and hole currents from fluxes
jn = jn_l + deltajn;
jp = jp_l + deltajp;

Jn = -jn*1000*par.e;
Jp = jp*1000*par.e;

% Total current
Jtot = Jn + Jp;


%Figures
if par.figson == 1
    
    % Open circuit voltage
    if par.OC == 1
        
        figure(7);
        plot (t, Voc);
        xlabel('Time [s]');
        ylabel('Voltage [V]');
        
    end
    
    % Dodgy way to change all the graphing but works!
    if par.OC == 1
        
        xnm(end) = round(xnm(end)/2);
        
    else
        
        xnm(end) = xnm(end);
    end
    %}
    %%%%% FIGURES %%%%%
    % Plotting defaults
    set(0,'DefaultLineLinewidth',1);
    set(0,'DefaultAxesFontSize',16);
    set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
    set(0,'DefaultAxesXcolor', [0, 0, 0]);
    set(0,'DefaultAxesYcolor', [0, 0, 0]);
    set(0,'DefaultAxesZcolor', [0, 0, 0]);
    set(0,'DefaultTextColor', [0, 0, 0]);
    
    for i=1:length(tarr)
        
        if pointtype == 't'
            timepoint = find(t <= tarr(i));
            pparr(i) = timepoint(end);
        elseif pointtype == 'V'
            Vpoint = find(solstruct.Vapp <= tarr(i));
            pparr(i) = Vpoint(end);
        end
        
        % Field
        Field = -gradient(V(pparr(i),:), x);
        
        % Calculates current at every point and all times -
        % UNRELIABLE FOR TOTAL CURRENT
        %if par.calcJ == 1
            
            [nloc,dnlocdx] = pdeval(0,x,n(pparr(i),:),x);
            [ploc,dplocdx] = pdeval(0,x,p(pparr(i),:),x);
            [iloc,dilocdx] = pdeval(0,x,a(pparr(i),:),x);
            [Vloc, dVdx] = pdeval(0,x,V(pparr(i),:),x);
            
            % Particle currents
            Jndiff = par.dev.mue.*par.kB*par.T.*dnlocdx*par.e;
            Jndrift = -par.dev.mue.*nloc.*dVdx*par.e;
            
            Jpdiff = -par.dev.muh*par.kB*par.T.*dplocdx*par.e;
            Jpdrift = -par.dev.muh.*ploc.*dVdx*par.e;
            
            Jidiff = -par.dev.muion*par.kB*par.T.*dilocdx*par.e;
            Jidrift = -par.dev.muion.*iloc.*dVdx*par.e;
            Jitot = Jidiff + Jidrift;
            
            % Particle current
            Jpart = Jndiff + Jndrift + Jpdiff + Jpdrift + Jidiff + Jidrift;
        
        
        % Displacement Current at right hand side
        Fend = -(dVdx(:, end));
        Jdispr = (par.e)*par.epp(3)*-gradient(dVdx(:, end), t);
        Jdispr = Jdispr';
        
        % TiO2 bulk
        % xrange = [200, 720];
        % PEDOT bulk
         xrange = [50, 570]'
        
        yrange = [-5e-7, 5e-7];
        % ion currents
        figure(20)
        plot(xnm, Jidiff)
        xlabel('Position [nm]');
        xlim([xrange(1), xrange(2)]);
        ylim([yrange(1), yrange(2)])
        ylabel('Ionic diffusion current [A]')
        hold on
        
        figure(21)
        plot(xnm, Jidrift)
        xlabel('Position [nm]');
        ylabel('Ionic drift current [A]')
        xlim([xrange(1), xrange(2)])
        ylim([yrange(1), yrange(2)])
        hold on

        figure(22)
        plot(xnm, Jitot)
        xlabel('Position [nm]');
        ylabel('Total ion current [A]')
        xlim([xrange(1), xrange(2)])
        ylim([yrange(1), yrange(2)])
        hold on
        
        figure(23)
        plot(xnm, (V(pparr(i), :)))
        xlabel('Position [nm]');
        ylabel('Electric field potential [V]')
        xlim([0, xnm(end)]);
        ylim([-0.5, 1.2])
        hold on
        
        figure(24)
        plot(xnm, Field)
        xlabel('Position [nm]');
        ylabel('Electric field [Vcm-1]')
        xlim([xrange(1), xrange(2)]);
        hold on
        %}
        
        %     % Pl intensity at time point
        %     PLint(pparr(i));
        %
        % Band Diagram
        FH1 = figure(1);
        %set(FigHandle, 'units','normalized','position',[.1 .1 .4 .4]);
        PH1 = subplot(3,1,1);
        plot (xnm, Efn(pparr(i),:), '--', xnm, Efp(pparr(i),:), '--', xnm, Ecb(pparr(i), :), xnm, Evb(pparr(i) ,:));
        %legend('E_{fn}', 'E_{fp}', 'CB', 'VB');
        set(legend,'FontSize',12);
        %xlabel('Position [nm]');
        ylabel('Energy [eV]');
        xlim([0, xnm(end)]);
        %ylim([-inf, 0.5]);
        set(legend,'FontSize',12);
        set(legend,'EdgeColor',[1 1 1]);
        grid off;
        drawnow;
        
        hold on
        
        % Final Charge Densities
        %figure(2)
        PH2 = subplot(3,1,2);
        semilogy(xnm, n(pparr(i), :), xnm, p(pparr(i), :));
        ylabel('{\itn, p} [cm^{-3}]')
        %legend('\itn', '\itp')
        %xlabel('Position [nm]')
        xlim([0, xnm(end)]);
        ylim([1e0, 1e20]);
        set(legend,'FontSize',12);
        set(legend,'EdgeColor',[1 1 1]);
        grid off
        
        hold on
        
        PH3 = subplot(3,1,3);
        plot(xnm, (rhoa(pparr(i),:))/1e18, 'black');
        ylabel('{\it\rho a} [x10^{18} cm^{-3}]');
        xlabel('Position [nm]');
        xlim([0, xnm(end)]);
        %ylim([0, 1.1*(max(sol(pparr(i),:,3))/1e19)]);
        set(legend,'FontSize',12);
        set(legend,'EdgeColor',[1 1 1]);
        grid off
        
        % %PL
        % figure(3)
        % plot(xnm, PL(pparr(i),:))
        % xlabel('Position [nm]');
        % ylabel('PL Intensity');
        
        % % ion plots
        % figure(3)
        % plot(xnm, a(pparr(i),:), xnm, Nionmat(pparr(i), :), xnm, rhoa(pparr(i), :))
        % xlabel('Position [nm]');
        % xlim([0, xnm(end)]);
        % legend('a', 'static', 'rho_a');
        
        % Current vs position
        figure(4)
        plot(xnm, Jn(pparr(i), :), xnm, Jp(pparr(i), :), xnm, Jtot(pparr(i), :))
        xlabel('Position [nm]')
        ylabel('Current density [mAcm-2]')
        hold on
        
        % figure(5)
        % plot(xnm, rhoc(end, :))
        % ylabel('Net Charge Density [cm^{-3}]')
        % xlabel('Position [nm]')
        % xlim([0, xnm(end)]);
        % set(legend,'FontSize',14);
        % set(legend,'EdgeColor',[1 1 1]);
        % grid off
        %
        % figure(6)
        % plot(t, ntot, t, ptot)
        % ylabel('Charge Density [cm^{-3}]')
        % xlabel('time [s]')
        % legend('electrons', 'holes')
        % set(legend,'FontSize',14);
        % set(legend,'EdgeColor',[1 1 1]);
        % grid off
        
        %
        if par.OM == 1 && par.Int~=0 || par.OM == 2 && par.Int~=0
            
            figure(7);
            plot(xnm, gx.AM15)
            ylabel('Generation Rate [cm^{3}s^{-1}]');
            xlabel('Position [nm]');
            legend('1 Sun');
            xlim([0, xnm(end)]);
            grid off
            
        end
        
        if par.calcJ == 1
            
            figure(8);
            plot(xnm,Jndiff(pparr(i), :),xnm,Jndrift(pparr(i), :),xnm,Jpdiff(pparr(i), :),xnm,Jpdrift(pparr(i), :),xnm,Jidiff(pparr(i), :),xnm,Jidrift(pparr(i), :),xnm,Jpart(pparr(i), :));
            legend('Jn diff','Jn drift','Jp diff','Jp drift','Ji diff','Ji drift','Total J');
            xlabel('Position [nm]');
            ylabel('Current Density [mA cm^-2]');
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            xlim([0, xnm(end)]);
            grid off;
            drawnow;
            
            hold on
        end
            %{
% Electric Field
figure(9);
surf(xnm, t, Floct);
xlabel('Position [m]');
ylabel('time [s]');
title('Electric Field');
            %}
            
    end
    
    figure(1)
    subplot(3,1,1);
    hold off
    subplot(3,1,2);
    hold off
    subplot(3,1,3);
    hold off
    figure(3)
    hold off
    figure(4)
    hold off
    figure(20)
    hold off
    figure(21)
    hold off    
    figure(22)
    hold off  
    figure(23)
    hold off
    figure(24)
    hold off
    
    
    if par.calcJ == 0 || par.calcJ == 1
        
        if par.JV == 1
            %JV
            figure(11)
            plot(Vapp_arr, Jtot(:, end))
            xlabel('V_{app} [V]')
            ylabel('Current Density [mA cm^-2]');
            grid off;
            
        else
            % Particle and displacement currents as a function of time
            figure(10);
            plot(t, Jtot(:, end));
            legend('Jtotal')
            xlabel('time [s]');
            ylabel('J [mA cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
            grid off;
            drawnow;
            
        end
        
        
    end
    
end

end


