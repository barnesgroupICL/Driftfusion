function calct(kscan, sol)

par = sol.p;
% device built-in
Vbi = par.PhiC - par.PhiA;
% kscan = scan rate [Vs-1]
Vapp = 0:0.001:1.4;
t = Vapp/kscan;
%t = 0.1:0.1:30;         % time [s]
Vapp = kscan.*t;    % Applied potential for plotting
alpha = 0.8;                % fraction of the potential across active layer;

deltat0 = t(2)-t(1);   % First time step to obtain F0
d = 510e-7;
F0 = alpha*kscan*deltat0/d;

Fint(1) = F0;

for i = 1:length(Vapp)
    % depletion width in contacts
    w1(i) = ((par.e*par.epp0*par.epp(1).*abs(-Vapp(i) + Vbi))/(par.e*par.NA(1))).^0.5;
    % limit depletion width at layer thickness
    if w1 >par.d(1)
        w1(i) = par.d(1);
    end
    w2(i) = ((par.e*par.epp0*par.epp(3).*abs(-Vapp(i) + Vbi))/(par.e*par.ND(3))).^0.5;
    if w2 >par.d(3)
        w2(i) = par.d(3);
    end
end
% Capacitance of junctions as a function of time:
C1 = (par.e*par.epp0*par.epp(1))./w1;
C2 = (par.e*par.epp0*par.epp(3))./w2;

Ctot = 1./(1./C1 + 1./C2);

%ionic resistance
Rion = d/(par.e*par.muion(2)*par.Nion(2));

tau = Ctot.*Rion;

for i = 2:length(t)-1
    deltat = t(i)-t(i-1);
    Fint(i) = (Fint(i-1)+alpha*kscan*deltat/d)*exp(-deltat./tau(i));
end

Finf = ((alpha*kscan*deltat/d)*exp(-deltat./tau))/(1-exp(-deltat./tau))

figure(200)
semilogy(Vapp(2:end), Fint)
xlabel('Vapp [V]')
ylabel('E [Vcm-1]')

% for i=1:length(V_thresh)
%
%     for j = 1:length(V_int0)
%
%         ratio_t1t2(i,j) = 3.*(log(V_thresh(i))-log(3*V_int0(j)))/(log(V_thresh(i))-log(V_int0(j)));
%
%     end
% end

end