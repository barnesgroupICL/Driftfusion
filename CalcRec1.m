function CalcRec1(struct, params)
% Calculate & plot recombination rates

set(0,'DefaultLineLinewidth',2);
set(0,'DefaultAxesFontSize',24);
set(0, 'DefaultFigurePosition', [600, 400, 640, 400]);

% Declare Variables - rememeber to add them here if adding to EParams
% The scoping rules for nested and anonymous functions require that all variables used within the function be present in the text of the code.

% Call parameters function
% VALUES MUST BE CORRECT FOR THE CHOSEN SOLUTION
%params = pnParams; 

% Unpacks params structure for use in current workspace 
v2struct(pnParams);  

sol = struct.sol;
t = struct.t*1e3;
x = struct.x;
xnm = x*1e7;

xpoints = length(x);

u1 = sol(:,:,1);
u2 = sol(:,:,2);

% p-type binary matrix
pBM = ones(length(t), xpoints)*diag(x <= tp | x>= xmax - tp);

pBA = ones(1, xpoints)*diag(x <= tp);% | x>= xmax - tp);
% Intrinsic binary matrix
%iBM = ones(length(t), xpoints)*diag(x > tp & x < (tp + ti) | x > (tp + ti + (2*tn)) & x <  xmax - tp);
% n-type binary matrix
nBM = ones(length(t), xpoints)*diag(x > tp & x <= (tp + 2*tn));

nBA = ones(1, xpoints)*diag(x > tp); %& x <= (tp + 2*tn));

%Usrh = zeros(length(t), length(x));
% Recomination Rate

for (i=1:xpoints/2)

Usrhm(i) = (kradhtl*max((u1(end, i)-htln0), (u2(end,i)-htlp0))).*pBA(i) + (kradetl*max((u1(end,i)-etln0), (u2(end,i)-etlp0))).*nBA(i);

end


%Usrhm = ((((u1.*u2)-ni^2)./((taun_htl*(u2+pthtl)) + (taup_htl*(u1+nthtl)))).*pBM) + ((((u1.*u2)-ni^2)./((taun_etl*(u2+ptetl)) + (taup_etl*(u1+ntetl)))).*nBM);

%Uradm = (kradetl*((u1.*u2)-(ni^2)).*pBM) + (krad*((u1.*u2)-(ni^2)).*iBM) + (kradhtl*((u1.*u2)-(ni^2)).*nBM);
%Uradm_act = (krad*((u1.*u2)-(ni^2)).*iBM);

%Utotm = (Usrhm + Uradm);

plot(xnm(1:floor(end/2)), Usrhm)
xlabel('Position [nm]');
ylabel('Rec rate [cm^{-3}]');
%{
Usrh = trapz(x, Usrhm, 2);                    % Divided by 2 for symmetric model     
Urad = trapz(x, Uradm, 2);
Urad_act = trapz(x, Uradm_act, 2);
Utot = trapz(x, Utotm, 2);
ntot = trapz(x, u1, 2);
ptot = trapz(x, u2, 2);

logUtot = log(Utot);
logntot = log(ntot);

% Reaction order
Delta = gradient(logUtot, logntot);

figure(50)
surf(x, t, Usrhm)
title ('SRH Recombination')
set(gca,'xlim',[0 (xmax/2)])

figure(51)
surf(x, t, Uradm)
title ('Radiative Recombination')
set(gca,'xlim',[0 (xmax/2)])

figure(52)
surf(x, t, Utotm)
title ('Total Recombination')
set(gca,'xlim',[0 (xmax/2)])

figure(53)
semilogy(t, Usrh, t, Urad, t, Utot)
legend('SRH', 'Rad Total', 'Total')
ylabel('Recombination Rate [cm^{-2}s^{-1}]');
xlabel('Time [s]');
set(legend,'FontSize',16);
set(legend,'EdgeColor',[1 1 1]);

figure(54)
semilogy(xnm, Uradm(1, :), xnm, Usrhm(1, :), xnm, Utotm(1, :), xnm, Uradm(end, :), '--', xnm, Usrhm(end, :), '--', xnm, Utotm(end, :), '--')
legend('Rad,i', 'SRH,i', 'Total,i','Rad,f', 'SRH,f', 'Total,f')
ylabel('Recombination Rate [cm^{-3}s^{-1}]');
xlabel('Position [nm]');
set(legend,'FontSize',14);
set(legend,'EdgeColor',[1 1 1]);
xlim([0, xmax*1e7/2]);

figure(56)
semilogy(t, ntot, t, ptot)
ylabel('n, p');
xlabel('Time [s]');
%}
rec_norm = Usrhm/max(Usrhm);

assignin('base', 'rec', Usrhm)
assignin('base', 'rec_norm', rec_norm)

end

