function F = makemovie(sol)
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultFigurePosition', [600, 400, 450, 300]);
set(0,'DefaultAxesXcolor', [0, 0, 0]);
set(0,'DefaultAxesYcolor', [0, 0, 0]);
set(0,'DefaultAxesZcolor', [0, 0, 0]);
set(0,'DefaultTextColor', [0, 0, 0]);
set(0,'DefaultLineLinewidth',1.5);

P = sol.p;

x = sol.x;
t = sol.t;
n = sol.sol(:,:,1);
p = sol.sol(:,:,2);
a = sol.sol(:,:,3);
V = sol.sol(:,:,4);

xnm = sol.x*1e7;


%% Binary matrices defining regions of the device

    % p-type binary matrix
    pBM = ones(length(t), P.xpoints)*diag(x <= P.dcum(1));
    % p-i interface binary matrix
    piBM = ones(length(t), P.xpoints)*diag(x > (P.dcum(1)) & x <= (P.dcum(1) + P.dint));       
    % Intrinsic binary matrix
    iBM = ones(length(t), P.xpoints)*diag(x > P.dcum(1) & x <= P.dcum(2));
    % i-n interface binary matrix II
    inBM = ones(length(t), P.xpoints)*diag(x > P.dcum(2) - P.dint & x <= P.dcum(2));
    % n-type binary matrix
    nBM = ones(length(t), P.xpoints)*diag(x > P.dcum(2) & x <= P.dcum(end));

    
% For graded bands
    % p-type binary matrix
    pBM2 = ones(length(t), P.xpoints)*diag(x <= P.dcum(1)-P.dint);
    % p-i interface binary matrix
    piBM2 = ones(length(t), P.xpoints)*diag(x > P.dcum(1) -P.dint & x <= P.dcum(1) + P.dint);       
    % Intrinsic binary matrix
    iBM2 = ones(length(t), P.xpoints)*diag(x > P.dcum(1) + P.dint & x <= P.dcum(2) - P.dint);
    % i-n interface binary matrix II
    inBM2 = ones(length(t), P.xpoints)*diag(x > P.dcum(2) - P.dint & x <= P.dcum(2) + P.dint);
    % n-type binary matrix
    nBM2 = ones(length(t), P.xpoints)*diag(x > P.dcum(2)  + P.dint & x <= P.dcum(end));

    
for k=1:P.tpoints

    j = 1;
    jj= 1;

    for i=1:P.xpoints
        
%         if P.OC == 1 && x(i) >= ceil(P.dcum(end)/2)
%             
%             i = P.xpoints - i + 1;
%             
%         end
        
        if x(i) > P.dcum(1) - P.dint && x(i) <= P.dcum(1)

            EApi1(k, i) = P.EA(1) + j*(P.dint/P.pint)*P.dEAdx(1);
            IPpi1(k, i) = P.IP(1) + j*(P.dint/P.pint)*P.dIPdx(1);
            N0pi1(k,i) = P.N0(1) + j*(P.dint/P.pint)*P.dN0dx(1);
            NIpi1(k,i) = 0 + j*(P.dint/P.pint)*(P.NI/(2*P.dint));
            j = j+1;

        else

            EApi1(k, i) = 0;
            IPpi1(k, i) = 0;
            N0pi1(k,i) = 0;
            NIpi1(k,i) = 0;
            
        end


        if x(i) > P.dcum(1) && x(i) <= P.dcum(1)+P.dint

            EApi2(k, i) = P.EA(1) + j*(P.dint/P.pint)*P.dEAdx(1);
            IPpi2(k, i) = P.IP(1) + j*(P.dint/P.pint)*P.dIPdx(1);
            N0pi2(k,i) = P.N0(1) + j*(P.dint/P.pint)*P.dN0dx(1);
            NIpi2(k,i) = 0 + j*(P.dint/P.pint)*(P.NI/(2*P.dint));
            j = j+1;

        else

            EApi2(k, i) = 0;
            IPpi2(k, i) = 0;
            N0pi2(k,i) = 0;
            NIpi2(k,i) = 0;
        end          
        
        
        if x(i) > P.dcum(2) -P.dint && x(i) <= P.dcum(2)

            EAin1(k, i) = P.EA(2) +jj*(P.dint/P.pint)*P.dEAdx(2);
            IPin1(k, i) = P.IP(2) +jj*(P.dint/P.pint)*P.dIPdx(2);
            N0in1(k,i) = P.N0(2) + jj*(P.dint/P.pint)*P.dN0dx(2);
            NIin1(k,i) = P.NI + jj*(P.dint/P.pint)*(-P.NI/(2*P.dint));
            jj = jj+1;

        else 

            EAin1(k, i) = 0;
            IPin1(k, i) = 0;
            N0in1(k,i) = 0;
            NIin1(k,i) = 0;
        end
        
        if x(i) > P.dcum(2) && x(i) <= P.dcum(2) + P.dint

            EAin2(k, i) = P.EA(2)+jj*(P.dint/P.pint)*P.dEAdx(2);
            IPin2(k, i) = P.IP(2)+jj*(P.dint/P.pint)*P.dIPdx(2);
            N0in2(k,i) = P.N0(2) + jj*(P.dint/P.pint)*P.dN0dx(2);
            NIin2(k,i) = P.NI + jj*(P.dint/P.pint)*(-P.NI/(2*P.dint));
            jj = jj+1;

        else %x(i) > P.dcum(2) +P.dint && x(i) <= P.dcum(2) + P.dint;

            EAin2(k, i) = 0;
            IPin2(k, i) = 0;
            N0in2(k,i) = 0;
            NIin2(k,i) = 0;
       end
        
    end
    
end  
    
%nstat = zeros(1, xpoints);                                  % Static charge array
nstat = (-P.NA(1)+P.ND(1))*pBM  + (-P.NA(2) + P.ND(2))*nBM + (-P.NA(3) + P.ND(3))*nBM;
rhoc = (-n + p + nstat);     % Net charge density calculated from adding individual charge densities

EA = P.EA(1)*pBM2 +  EApi1 + EApi2 + P.EA(2)*iBM2 + EAin1 + EAin2 + P.EA(3)*nBM2;
IP = P.IP(1)*pBM2 + IPpi1 + IPpi2 + P.IP(2)*iBM2 + IPin1 + IPin2 + P.IP(3)*nBM2;
N0 = P.N0(1)*pBM2  + N0pi1 + N0pi2 + P.N0(2)*iBM2 + N0in1 + N0in2 +P.N0(3)*nBM2;
NImat = NIpi1 + NIpi2 + P.NI*iBM2 + NIin1 + NIin2;
Ei = P.Eif(1)*pBM  + P.Eif(2)*iBM + P.Eif(3)*nBM;
ni = P.ni(1)*pBM  + P.ni(2)*iBM + P.ni(3)*nBM;

Ecb = EA-V;                                 % Conduction band potential
Evb = IP-V;                                 % Valence band potential
Efn = real(Ecb+(P.kB*P.T/P.q)*log(n./N0));        % Electron quasi-Fermi level 
Efp = real(Evb-(P.kB*P.T/P.q)*log(p./N0));        % Hole quasi-Fermi level
Phin = real(Ei+(P.kB*P.T/P.q)*log(n./ni)-EA);     % Chemical Potential electrons
Phip = real(Ei-(P.kB*P.T/P.q)*log(p./ni)-EA);     % Chemical Potential holes
Phi = Phin - Phip;

for i=1:length(sol.t)

fig1 = figure(30)
clf
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
Vnow = round(sol.Vapp(i), 2, 'decimal');
anno = ['V = ', num2str(Vnow), ' V'];
T = annotation('textbox', dim, 'String', anno, 'FitBoxToText','on');
T.FontSize = 16;

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
%{
PH2 = subplot(3,1,3);
plot(xnm, a(i, :))
xlabel('Position')
ylabel('Ion density [cm^{-3}]')
xlim([1, xnm(end)])
%ylim([1e6, 1e18])
%}

% semilogy(xnm, Usrh(i, :))
% xlabel('Position [nm]')
% ylabel('Urec [cm^{-3}s^{-1}]')
% xlim([1, xnm(end)])
% ylim([1e21, 1e25])

F(i) = getframe(fig1);

%
end


end
