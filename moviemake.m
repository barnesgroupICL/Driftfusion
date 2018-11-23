function F = moviemake(solstruct)

% Makes a frame file F from a solution
% Currently configured to output the band diagram

% Plotting defaults
set(0,'DefaultLineLinewidth',1);
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
set(0,'DefaultAxesXcolor', [0, 0, 0]);
set(0,'DefaultAxesYcolor', [0, 0, 0]);
set(0,'DefaultAxesZcolor', [0, 0, 0]);
set(0,'DefaultTextColor', [0, 0, 0]);

ionfigon = 0;
capfigon = 0;

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
taunmat = repmat(par.dev.taun, length(t), 1);
taupmat = repmat(par.dev.taup, length(t), 1);
ntmat =  repmat(par.dev.nt, length(t), 1);
ptmat =  repmat(par.dev.pt, length(t), 1);

Ecb = EAmat-V;                                 % Conduction band potential
Evb = IPmat-V;                                 % Valence band potential

Efn = zeros(size(n,1), size(n,2));
Efp = zeros(size(n,1), size(n,2));

% charge density
rho = -n + p + a -NAmat + NDmat - Nionmat;

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
    dadt(:,j) = gradient(a(:,j), t);
end

dndtInt = trapz(x, dndt, 2);
dpdtInt = trapz(x, dpdt, 2);
dadtInt = trapz(x, dadt, 2);
% Recombination
Ubtb = kradmat.*(n.*p - nimat.^2);
% Recombination
Ubtb = kradmat.*(n.*p - nimat.^2);

Usrh = ((n.*p - nimat.^2)./((taunmat.*(p+ptmat)) + (taupmat.*(n+ntmat))));

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
djadx = -dadt;

% Integrate across the device to get delta fluxes at all positions
deltajn = cumtrapz(par.x, djndx, 2);
deltajp = cumtrapz(par.x, djpdx, 2);
deltaja = cumtrapz(par.x, djadx, 2);
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
ja = 0 + deltaja;

Jn = -jn*par.e;
Jp = jp*par.e;
Ja = ja*par.e;

% Total electronic current
Jtot = Jn + Jp;

j = 1;
k = 1;
% Iteration around number of layers and desired points
for i=1:2*length(par.parr)-1
    
    if rem(i, 2) == 1
        parrint(j) = par.parr(k);
        j = j+1;
        k = k+1;
    elseif rem(i, 2) == 0
        parrint(j) = par.pint;
        j = j+1;
    end
end
pcum = cumsum(parrint);
pcum = [0,pcum]+1;
pcum(end) = pcum(end)-1;

% TiO2 bulk
xrange = [0, xnm(end)];

% xrange = [200, 724];
% PEDOT bulk
% xrange = [50, 574]'
% Tio2, Spiro-pvsk interface
% xrange = [198, 210];
% PEDOT-pvsk interface
% xrange = [48, 60];
pp1 = find(Vapp_arr < 1);
pp1 = pp1(end);

for i=1:pp1%length(Vapp_arr)
    
    fig1 = figure(30)
    clf

    % Field
    Field = -gradient(V(i,:), x);
    
    % Calculates current at every point and all times -
    % UNRELIABLE FOR TOTAL CURRENT
    %if par.calcJ == 1
    
    [nloc,dnlocdx] = pdeval(0,x,n(i,:),x);
    [ploc,dplocdx] = pdeval(0,x,p(i,:),x);
    [iloc,dilocdx] = pdeval(0,x,a(i,:),x);
    [Vloc, dVdx] = pdeval(0,x,V(i,:),x);
    
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
    %{
    yrange1 = [-2e-7, 1e-7];
    yrange2 = [-1e-7, 2e-7];
    
    subplot(2,1,1);
    plot(xnm, Jidiff)
    xlabel('Position [nm]');
    xlim([xrange(1), xrange(2)]);
    ylim([yrange1(1), yrange1(2)])
    ylabel('Ionic diffusion current [A]')
    
    subplot(2,1,2);
    plot(xnm, Jidrift)
    xlabel('Position [nm]');
    ylabel('Ionic drift current [A]')
    xlim([xrange(1), xrange(2)])
    ylim([yrange2(1), yrange2(2)])
    
    dim = [.15 0.15 .3 .3];
    tnow = round(solstruct.t(i), 2, 'decimal');
    anno = ['t = ', num2str(tnow), ' s'];
    T = annotation('textbox', dim, 'String', anno, 'FitBoxToText','on');
    T.FontSize = 16;
    %}
PH1 = subplot(3,1,1);
plot (xnm, Efn(i,:), '--', xnm, Efp(i,:), '--', xnm, Ecb(i, :), xnm, Evb(i ,:));
%legend('E_{fn}', 'E_{fp}', 'CB', 'VB');
set(legend,'FontSize',12);
xlabel('Position [nm]');
ylabel('Energy [eV]');
xlim([0, xnm(end)]);
ylim([-7, -3]);
set(legend,'FontSize',12);
set(legend,'EdgeColor',[1 1 1]);
grid off;
drawnow;

dim = [.2 0 .3 .3];
Vnow = round(solstruct.Vapp(i), 2, 'decimal');
%tnow = round(solstruct.t(i), 2, 'decimal');
anno = ['V = ', num2str(Vnow), ' V'];
%anno = ['t = ', num2str(tnow), ' s'];
T = annotation('textbox', dim, 'String', anno, 'FitBoxToText','on');
T.FontSize = 16;

PH2 = subplot(3,1,2);
semilogy(xnm, n(i, :), xnm, p(i,:))
xlabel('Position [nm]')
ylabel('Carrier density [cm^{-3}]')
xlim([0, xnm(end)])
ylim([1e6, 1e19])

%PH1 = subplot(3,1,1);
% plot(xnm, -V(i, :))
% xlim([1, xnm(end)])
% ylim([-1.4, 0])
% xlabel('Position [nm]')
% ylabel('Electric potential [V]')

%dim = [.2 .5 .3 .3];
%anno = ['V = ', num2str(sol.Vapp(i))];
%T = annotation('textbox', dim, 'String', anno, 'FitBoxToText','on');

%PH2 = subplot(3,1,2);
% semilogy(xnm, n(i, :), xnm, p(i,:))
% xlabel('Position [nm]')
% ylabel('Carrier density [cm^{-3}]')
% xlim([1, xnm(end)])
% ylim([1e6, 1e18])

PH2 = subplot(3,1,3);
plot(xnm, rhoa(i, :))
xlabel('Position')
ylabel('Ion density [cm^{-3}]')
xlim([0, xnm(end)])
ylim([-1e18, 4e18])

%{
        % electron and hole currents as function of position from continuity
        %figure(110)
        plot(xnm, Jn(i, :), xnm, Jp(i, :), xnm, Jtot(i, :), '--')
        xlabel('Position [nm]')
        ylabel('J [A]')    
        xlim([xrange(1), xrange(2)])
        ylim([-8, 8])
        %legend('Jn', 'Jp', 'Jtot')
        %hold on
%}
% semilogy(xnm, Usrh(i, :))
% xlabel('Position [nm]')
% ylabel('Urec [cm^{-3}s^{-1}]')
% xlim([1, xnm(end)])
% ylim([1e21, 1e25])
    %}
    F(i) = getframe(fig1);
    
    %
end


end
