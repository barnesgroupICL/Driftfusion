function [Vapp_arr, Jtot] = pinana(varargin)

% Plotting defaults
set(0,'DefaultLineLinewidth',1);
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
set(0,'DefaultAxesXcolor', [0, 0, 0]);
set(0,'DefaultAxesYcolor', [0, 0, 0]);
set(0,'DefaultAxesZcolor', [0, 0, 0]);
set(0,'DefaultTextColor', [0, 0, 0]);

ionfigon = 0;
ioncur = 0;
capfigon = 0;
cardent = 0;
recfigon = 0;
ddcur = 0;
transfigon =1;

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
    
    Vapp_arr = Efn(:, end) - Efp(:, 1);
    
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

Usrh = ((n.*p - nimat.^2)./((taunmat.*(p+ptmat)) + (taupmat.*(n+ntmat))));

U = Ubtb + Usrh;

Usrhnogen = ((n.*p)./((taunmat.*(p+ptmat)) + (taupmat.*(n+ntmat))));

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

% Total current
Jtot = Jn + Jp + Ja;

%% Capacitance of the interfaces
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
for i=1:length(t)
    % Charge densities across interfaces
    rhoint1(i) = trapz(x(pcum(1):pcum(2)-1), rho(i, pcum(1):pcum(2)-1))*par.e;
    rhoion1(i) = trapz(x(pcum(2):(pcum(3)-1)+round(parrint(3)/2)), rho(i, pcum(2):(pcum(3)-1)+round(parrint(3)/2)))*par.e;
    rhoint2(i) = trapz(x(pcum(5):pcum(6)-1), rho(i, pcum(5):pcum(6)-1))*par.e;
    rhoion2(i) = trapz(x(pcum(3)+round(parrint(3)/2):pcum(5)-1), rho(i, pcum(3)+round(parrint(3)/2):pcum(5)-1))*par.e;
end

% Voltage drop across the interfaces
Vint1 = -(V(:, pcum(1)) - V(:, pcum(3)+30));
Vint2 = -(V(:, pcum(4)-30) - V(:, end));

% Voltage drop within the perovskite
Vint1_per = -(V(:, pcum(2)) - V(:, pcum(3)+30));
Vint2_per = -(V(:, pcum(4)-30) - V(:, pcum(5)));

if ioncur == 1
    %% ion current at centre of active layer as a function of time
    Ja_mid = Ja(:, pcum(3)-1+round(parrint(3)/2));
    logJa_mid = log10(Ja_mid/max(Ja_mid));
    
    p1 = find(logJa_mid < -0.1);
    p1 = p1(1);
    
    p2 = find(logJa_mid < -4);
    p2 = p2(1);
    
    fitob = fit(t(p1:p2)', logJa_mid(p1:p2), 'poly1');
    
    % normalised ion current in centre of active layer as function of time
    figure(2)
    plot(fitob, t(p1:p2), logJa_mid(p1:p2))
    xlabel('Time [s]')
    ylabel('log10(Ja) [norm]')
    
    fitob
    
    figure(40)
    plot(t, Ja_mid)
    xlabel('Time [s]')
    ylabel('Ja [A]')
    
end

% Electric field
for i = 1:length(t)
    Field(i, :) = gradient(V(i, :), x);
end

% Vrec - this doesn't really work as we need the maxima and minima
% TiO2
Vrec1_0 = Ecb(1, pcum(3))- Efn(1, pcum(3));
Vrec2_0 = Efp(1, pcum(4)-1) - Evb(1, pcum(4)-1);
VRec0_verf = Vrec1_0 + Vrec2_0;

Vrec1 = Vrec1_0-(Ecb(:, pcum(3))- Efn(:, pcum(3)));
Vrec2 = Vrec2_0-(Efp(:, pcum(4)-1) - Evb(:, pcum(4)-1)); 
Efint = Efn(:, pcum(3)) - Efp(:, pcum(4)-1);
Everf = Vrec1 + Vrec2;
Jrec1 = Jn(:, pcum(3));
Jrec2 = Jp(:, pcum(4)-1);

% Vinj
Vinj1 = (Efp(1, pcum(3))-Evb(1, pcum(3)))-(Efp(:, pcum(3))-Evb(:, pcum(3)));
Vinj2 = (Ecb(1, pcum(4)-1)- Efn(1, pcum(4)-1))-(Ecb(:, pcum(4)-1)- Efn(:, pcum(4)-1));
Jinj1 = Jp(:, pcum(3));
Jinj2 = Jn(:, pcum(4));

Jbulk = trapz(x(pcum(3):pcum(4)-1), U(pcum(3):pcum(4)-1));


%Figures
if par.figson == 1
    
    rs = 0;
    %% Specify ranges for plotting- should be an input argument
    % whole device
    xrange = [xnm(pcum(1)), xnm(pcum(end))];
    % bulk
    % xrange = [xnm(pcum(2)), xnm(pcum(5))];
    % Interface 1:
    %xrange = [xnm(pcum(3))-10, xnm(pcum(3))+6];
    %Interface 2:
    %xrange = [xnm(pcum(4))-6, xnm(pcum(4))+10];
    
    yrange = [-inf, inf];
    % yrange = [-5e-7, 5e-7];
    % ion currents
    
    %%%%% FIGURES %%%%%
    % Plotting defaults
    set(0,'DefaultLineLinewidth',1);
    set(0,'DefaultAxesFontSize',16);
    set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
    set(0,'DefaultAxesXcolor', [0, 0, 0]);
    set(0,'DefaultAxesYcolor', [0, 0, 0]);
    set(0,'DefaultAxesZcolor', [0, 0, 0]);
    set(0,'DefaultTextColor', [0, 0, 0]);
    
    
    % Field strength at centre of active layer and mean field
    % strength across perovskite as function of time
%     figure(2)
%     plot(Vapp_arr(2:end), -Field((2:end), pcum(3)+round(parrint(3)/2)));%, t, mean(Field(:, pcum(3):pcum(4)), 2))
%     xlabel('Vapp [V]')
%     ylabel('Electric field [Vcm-1]')
    %legstr = [par.stack(1), num2str(round(mean(Vapp_arr/t)), 2),'Vs-1'];
    %legend([par.stack(1), num2str(round(mean(Vapp_arr/t))),'Vs-1'])
    
    
    for i=1:length(tarr)
        
        if pointtype == 't'
            timepoint = find(t <= tarr(i));
            pparr(i) = timepoint(end);
        elseif pointtype == 'V'
            
            if rs ==1
                Vpoint = find(solstruct.Vapp <= tarr(i));
                pparr(i) = Vpoint(1);
            else
                Vpoint = find(solstruct.Vapp <= tarr(i));
                pparr(i) = Vpoint(end);
            end
        end
        
        % electron and hole currents as function of position from continuity
        figure(110)
        plot(xnm, Jn(pparr(i), :), xnm, Jp(pparr(i), :), xnm, Jtot(pparr(i), :))
        xlabel('Position [nm]')
        ylabel('J [A]')    
        xlim([xrange(1), xrange(2)])
        hold on

        if transfigon == 1
            
            figure(90)
            plot(Vrec1, Jrec1, Vrec2, Jrec2)
            xlabel('Vrec [V]')
            ylabel('Jrec [A]')
            legend('Jrec1', 'Jrec2')
            
            lnJ = log(Jtot(end));
            lnJrec1 = log(Jrec1);
            lnJrec2 = log(Jrec2);
            
            lnJinj1 = log(Jinj1);
            lnJinj2 = log(Jinj2);
            
            p1 = find(Vrec1 > 0.3);
            p1 = p1(1);
            
            p2 = find(Vrec1 > 0.4);
            p2 = p2(1);
            
            fitJrec1 = fit(Vrec1(p1:p2), lnJrec1(p1:p2), 'poly1');
            
            p1 = find(Vrec2 > 0.2);
            p1 = p1(1);
            
            p2 = find(Vrec2 > 0.3);
            p2 = p2(1);           
            
            fitJrec2 = fit(Vrec2(p1:p2), lnJrec2(p1:p2), 'poly1');
            
            % normalised ion current in centre of active layer as function of time
            figure(91)
            plot(Vrec1, lnJrec1)
            xlabel('Vrec1 [V]')
            ylabel('ln(Jrec)')
            legend('Jrec1')
%             dim = [.2 .5 .3 .3];
%             anno = ['gradient =' num2str(round(fitJrec1.p1, 2)), newline, 'intercept =', num2str(round(fitJrec1.p2, 2))];
%             annotation('textbox',dim,'String',anno,'FitBoxToText','on', 'Fontsize', 16);
            
            figure(92)
            plot(Vrec2, lnJrec2)
            xlabel('Vrec2 [V]')
            ylabel('ln(Jrec)')
            legend('Jrec2')
%             dim = [.2 .5 .3 .3];
%             anno = ['gradient =' num2str(round(fitJrec2.p1, 2)), newline, 'intercept =', num2str(round(fitJrec2.p2, 2))];
%             annotation('textbox',dim,'String',anno,'FitBoxToText','on', 'Fontsize', 16);
 
            figure(93)
            plot(Vapp_arr, lnJrec1, Vapp_arr, lnJrec2, Vapp_arr, lnJinj1, Vapp_arr, lnJinj2)
            xlabel('Vapp [V]')
            ylabel('ln(Jinj)')
            legend('Jrec1', 'Jrec2', 'Jinj1', 'Jinj2')
%             dim = [.2 .5 .3 .3];
%             anno = ['gradient =' num2str(round(fitJinj1.p1, 2)), newline, 'intercept =', num2str(round(fitJinj1.p2, 2))];
%             annotation('textbox',dim,'String',anno,'FitBoxToText','on', 'Fontsize', 16);
            
            figure(94)
            plot(Vinj2, lnJinj2)
            xlabel('Vinj2 [V]')
            ylabel('ln(Jinj)')
            legend('Jinj2')
%             dim = [.2 .5 .3 .3];
%             anno = ['gradient =' num2str(round(fitJinj2.p1, 2)), newline, 'intercept =', num2str(round(fitJinj2.p2, 2))];
%             annotation('textbox',dim,'String',anno,'FitBoxToText','on', 'Fontsize', 16);
            
            figure(95)
            plot(Vapp_arr, Vrec1, Vapp_arr, Vrec2)
            xlabel('Vapp [V]')
            ylabel('Vrec [V]')
            xlim([Vapp_arr(1), 1])%Vapp_arr(end)])
            legend('Vrec1', 'Vrec2')  
            
            figure(9)
            plot(Vapp_arr, log(Jtot(:, end)))
            xlabel('Vapp [V]')
            ylabel('ln(J)')
            xlim([Vapp_arr(1), Vapp_arr(end)])
            
        end
        
        
        if recfigon == 1
            
            figure(60)
            semilogy(xnm, Usrh(pparr(i), :));%, xnm, Ubtb(pparr(i),:), xnm, U(pparr(i),:))
            xlabel('Position [nm]');
            ylabel('Rec rate [cm-3s-1]')
            %legend('SRH', 'rad', 'tot')
            xlim([xrange(1), xrange(2)])
            ylim([yrange(1), yrange(2)])
            hold on
            
            for j=1:length(t)
                Usrh_bulk(j) = trapz(x(pcum(3):pcum(4)-1), Usrh(j, pcum(3):pcum(4)-1),2);
                Usrh_int1(j) = trapz(x(pcum(2):pcum(3)-1), Usrh(j, pcum(2):pcum(3)-1),2);
                Usrh_int2(j) = trapz(x(pcum(4):pcum(5)-1), Usrh(j, pcum(4):pcum(5)-1),2);
                
                Usrhnogen_int1(j) = trapz(x(pcum(2):pcum(3)-1), Usrhnogen(j, pcum(2):pcum(3)-1),2);
                Usrhnogen_int2(j) = trapz(x(pcum(4):pcum(5)-1), Usrhnogen(j, pcum(4):pcum(5)-1),2);
                
                Ubtb_bulk(j) = trapz(x(pcum(3):pcum(4)-1), Ubtb(j, pcum(3):pcum(4)-1),2);
                Ubtb_int1(j) = trapz(x(pcum(2):pcum(3)-1), Ubtb(j, pcum(2):pcum(3)-1),2);
                Ubtb_int2(j) = trapz(x(pcum(4):pcum(5)-1), Ubtb(j, pcum(4):pcum(5)-1),2);
            end
            % srh as a function off time
            
            if par.JV ==1
                figure(61)
                semilogy(Vapp_arr, Usrh_bulk, Vapp_arr, Usrh_int1, Vapp_arr, Usrh_int2)
                xlabel('Vapp [V]')
                ylabel('SRH rate [cm-2s-1]')
                legend('bulk', 'int1', 'int2')
                
                if rs == 0
                    xlim([Vapp_arr(1), Vapp_arr(end)])
                else
                    xlim([Vapp_arr(end), Vapp_arr(1)])
                end
                
                %                 % delta srh
                %                 figure(62)
                %                 semilogy(Vapp_arr, (Usrh_bulk-Usrh_bulk(1)), Vapp_arr, (Usrh_int1-Usrh_int1(1)), Vapp_arr, (Usrh_int2-Usrh_int2(1)))
                %                 xlabel('Vapp [V]')
                %                 ylabel('Delta SRH rate [cm-3s-1]')
                %                 legend('bulk', 'int1', 'int2')
                
                %                 if rs == 0
                %                     xlim([Vapp_arr(1), Vapp_arr(end)])
                %                 else
                %                     xlim([Vapp_arr(end), Vapp_arr(1)])
                %                 end
                
                
                
            else
                
                figure(61)
                plot(t, Usrh_int1*par.e, t, Ubtb_int1*par.e, t, (Usrh_int1+Ubtb_int1)*par.e)
                xlabel('Time [s]')
                ylabel('Interface 1 rec current [Acm-2s-1]')
                legend('SRH', 'btb', 'total')
                
                figure(62)
                plot(t, Usrh_int2*par.e, t, Ubtb_int2*par.e, t, (Usrh_int2+Ubtb_int2)*par.e)
                xlabel('Time [s]')
                ylabel('Interface 2 rec current [Acm-2s-1]')
                legend('SRH', 'btb', 'total')
                
                figure(63)
                plot(t, Usrhnogen_int1*par.e, t, Usrhnogen_int2*par.e)
                xlabel('Time [s]')
                ylabel('SRH rec current [Acm-2s-1]')
                legend('int1', 'int2', 'total')
                %                 % delta srh
                %                 figure(62)
                %                 semilogy(t, (Usrh_bulk-Usrh_bulk(1)), t, (Usrh_int1-Usrh_int1(1)), t, (Usrh_int2-Usrh_int2(1)))
                %                 xlabel('Time [s]')
                %                 ylabel('Delta SRH rate [cm-3s-1]')
                %                 legend('bulk', 'int1', 'int2')
                %
                
            end
        end
        
        
        
        if capfigon == 1
            
            %             % charge density as a function of position
            %             figure(30)
            %             plot(xnm, rho(pparr(i), :))
            %             xlabel('Position [nm]');
            %             ylabel('Charge density [cm-3]')
            %             hold on
            %
            % capacitance as a function of voltage
            Cint1 = gradient(rhoint1, Vint1);
            Cion1 = gradient(rhoion1, Vint1);
            Cint1 = abs(Cint1);
            Cion1 = abs(Cion1);
            
            Cint2 = gradient(rhoint2, Vint2);
            Cion2 = gradient(rhoion2, Vint2);
            Cint2 = abs(Cint2);
            Cion2 = abs(Cion2);
            
            Cinttot = 1./(1./Cint1 + 1./Cint2);
            Ciontot = 1./(1./Cion1 + 1./Cion2);
            
            % ionic resistance
            a_av = mean(a(:, pcum(3):pcum(4)),2);
            sigma_ion = (par.e.*par.muion(2).*a_av);
            Rion = par.d(2)./sigma_ion;
            
            % RC time constants
            tau_RCint = Rion'.*Cinttot;
            tau_RCion = Rion'.*Ciontot;
            
            figure(31)
            semilogy(t, Cint1, t, Cion1)
            xlabel('Voltage [V]')
            ylabel('C = dQ/dV [Fcm-2]')
            ylim([1e-9, 1e-5])
            xlim([t(1), t(end)])
            legend('Cint1', 'Cion1')
            
            figure(32)
            semilogy(t, Cint2, t, Cion2)
            xlabel('Voltage [V]')
            ylabel('C = dQ/dV [Fcm-2]')
            ylim([1e-9, 1e-5])
            xlim([t(1), t(end)])
            legend('Cint2', 'Cion2')
            
            % Voltage drop across perovskite
            figure(33)
            plot(t, Vint1_per, t, Vint2_per)
            xlabel('Applied Voltage [V]')
            ylabel('Voltage drop in perovskite [V]')
            legend('Vint1', 'Vint2')
            
            % Voltage drop across interfaces
            figure(34)
            plot(t, Vint1, t, Vint2)
            xlabel('Applied Voltage [V]')
            ylabel('Voltage drop across interfaces [V]')
            legend('Vint1', 'Vint2')
            
            figure (35)
            plot(Vapp_arr, Ciontot)
            xlabel('Voltage [V]')
            ylabel('Total capacitance [Fcm-2]')
            %ylim([1e-9, 1e-5])
            xlim([Vapp_arr(1), Vapp_arr(end)])
            %legend('Cint', 'Cion')
            hold on
            figure(36)
            plot(Vapp_arr, tau_RCion)
            xlabel('Voltage [V]')
            ylabel('RC time constant [s]')
            xlim([Vapp_arr(1), Vapp_arr(end)])
            %legend('tau int', 'tau ion')
            hold on
            
        end
        
        
        % Field
        Fieldpoint = -gradient(V(pparr(i),:), x);
        
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
        
        
        if ddcur == 1
            
            figure(80)
            plot(xnm, Jndrift, xnm, Jndiff)
            xlabel('Position [nm]');
            xlim([xrange(1), xrange(2)]);
            ylim([yrange(1), yrange(2)])
            ylabel('Jn [A]')
            legend('Drift', 'Diff')
            hold on
            
            figure(81)
            plot(xnm, Jpdrift, xnm, Jpdiff)
            xlabel('Position [nm]');
            xlim([xrange(1), xrange(2)]);
            ylim([yrange(1), yrange(2)])
            ylabel('Jp [A]')
            legend('Drift', 'Diff')
            hold on
        end
        
        
        
        if ionfigon ==1
            
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
            
            %             figure(22)
            %             plot(xnm, Jitot)
            %             xlabel('Position [nm]');
            %             ylabel('Total ion current [A]')
            %             %xlim([xrange(1), xrange(2)])
            %             %ylim([yrange(1), yrange(2)])
            %             hold on
            
            figure(23)
            plot(xnm, (-V(pparr(i), :)))
            xlabel('Position [nm]');
            ylabel('-Electric field potential [V]')
            xlim([0, xnm(end)]);
            ylim([-1.2, 0.5])
            hold on
            
            figure(24)
            plot(xnm, Fieldpoint)
            xlabel('Position [nm]');
            ylabel('Electric field [Vcm-1]')
            xlim([xrange(1), xrange(2)]);
            hold on
            
            figure(3)
            plot(xnm, rhoa(pparr(i), :)) %xnm, a(pparr(i),:), xnm, Nionmat(pparr(i), :),
            ylabel('Ionic space charge density [cm-3]')
            xlabel('Position [nm]');
            xlim([xrange(1), xrange(2)]);
            %legend('a', 'static', 'rho_a');
            hold on
            
            % Potential across the active layer as a function of time
            p1 = find(xnm <= xrange(1));
            p1 = p1(end);
            p2 = find(xnm <= xrange(2));
            p2 = p2(end);
            
            figure(25)
            plot(t, -((V(:,p2)-V(:,p1))-(V(end,p2)-V(end,p1))), t, -((V(:,end)-V(:,1))-(V(end,end)-V(end,1))))
            legend('active layer', 'device')
            xlabel('time [s]')
            ylabel('Delta V [V]')
            
            % total ionic current
            figure(26)
            plot(xnm, Ja(pparr(i), :))
            xlabel('Position [nm]');
            ylabel('Total ion current [A]')
            xlim([xrange(1), xrange(2)]);
            %xlim([xnm(1), xnm(end)])
            %ylim([yrange(1), yrange(2)])
            hold on
            
        end
        
        % Carrier densities at interfaces as a function of time
        if cardent == 1
            
            figure(50)
            semilogy(Vapp_arr, n(:, pcum(3)), Vapp_arr, p(:, pcum(3)))
            xlabel('Time [s]')
            ylabel('n,p at HTL/pvsk interface')
            legend('n','p')
            
            figure(51)
            semilogy(Vapp_arr, n(:, pcum(4)), Vapp_arr, p(:, pcum(4)))
            xlabel('Time [s]')
            ylabel('n,p at pvsk/ETL interface')
            legend('n','p')
        end
        
        %     % Pl intensity at time point
        %     PLint(pparr(i));
        %
        % Band Diagram
        %FH1 = figure(1);
        %set(FigHandle, 'units','normalized','position',[.1 .1 .4 .4]);
        %PH1 = subplot(3,1,1);
        figure(70)
        plot (xnm, Efn(pparr(i),:), '--', xnm, Efp(pparr(i),:), '--', xnm, Ecb(pparr(i), :), xnm, Evb(pparr(i) ,:));
        %legend('E_{fn}', 'E_{fp}', 'CB', 'VB');
        set(legend,'FontSize',12);
        xlabel('Position [nm]');
        ylabel('Energy [eV]');
        xlim([xrange(1), xrange(2)]);
        %xlim([0, xnm(end)]);
        %ylim([-inf, 0.5]);
        set(legend,'FontSize',12);
        set(legend,'EdgeColor',[1 1 1]);
        grid off;
        drawnow;
        
        hold on
        
        % Final Charge Densities
        %figure(2)
        %PH2 = subplot(3,1,2);
        figure(71)
        semilogy(xnm, n(pparr(i), :), xnm, p(pparr(i), :));
        ylabel('{\itn, p} [cm^{-3}]')
        %legend('\itn', '\itp')
        xlabel('Position [nm]')
        %xlim([0, xnm(end)]);
        xlim([xrange(1), xrange(2)]);
        ylim([1e0, 1e20]);
        set(legend,'FontSize',12);
        set(legend,'EdgeColor',[1 1 1]);
        grid off
        
        hold on
        
        %PH3 = subplot(3,1,3);
        figure(72)
        plot(xnm, (rhoa(pparr(i),:))/1e18, 'black');
        ylabel('{\it\rho a} [x10^{18} cm^{-3}]');
        xlabel('Position [nm]');
        %xlim([0, xnm(end)]);
        xlim([xrange(1), xrange(2)]);
        %ylim([0, 1.1*(max(sol(pparr(i),:,3))/1e19)]);
        set(legend,'FontSize',12);
        set(legend,'EdgeColor',[1 1 1]);
        grid off
        
        % %PL
        % figure(3)
        % plot(xnm, PL(pparr(i),:))
        % xlabel('Position [nm]');
        % ylabel('PL Intensity');
        
        %         % Current vs position
        %         figure(4)
        %         plot(xnm, Jn(pparr(i), :), xnm, Jp(pparr(i), :), xnm, Jtot(pparr(i), :))
        %         xlabel('Position [nm]')
        %         ylabel('Current density [Acm-2]')
        %         hold on
        
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
            ylabel('Current Density [A cm^-2]');
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
    
    %     figure(1)
    %     subplot(3,1,1);
    %     hold off
    %     subplot(3,1,2);
    %     hold off
    %     subplot(3,1,3);
    %     hold off
    figure(70)
    hold off
    figure(71)
    hold off
    figure(72)
    hold off
    
    if ionfigon ==1
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
        figure(26)
        hold off
    end
    
    if recfigon == 1
        figure(60)
        %hold off
    end
    
    if capfigon == 1
        figure(30)
        hold off
    end
    
    
    if par.calcJ == 0 || par.calcJ == 1
        
        if par.JV == 1
            %JV
            figure(11)
            plot(Vapp_arr, Jtot(:, end))
            xlabel('V_{app} [V]')
            ylabel('Current Density [A cm^-2]');
            grid off;
            
        else
            % Particle and displacement currents as a function of time
            figure(10);
            plot(t, Jtot(:, end));
            legend('Jtotal')
            xlabel('time [s]');
            ylabel('J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
            grid off;
            drawnow;
            
        end
        
        
    end
    
end

end


