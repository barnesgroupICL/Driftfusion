function nidt = transient_nid_ana(sol_OC)
% Takes an input from transient_nid and plots nid(t)
par = sol_OC(1).par;
% Incident photon flux density at 1 Sun across device
% Get x_halfi
xihalf = getvarihalf(sol_OC(1).x);
G = trapz(xihalf, par.gx1);
   
% Extract desrired values from the solution
for i = 1:length(sol_OC)
   % Split solution
   [u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol_OC(i));
   
   Voct(i,:) = dfana.Voct(sol_OC(i));
   t = sol_OC(i).t;     
   Int(i) = sol_OC(i).par.int1;
   % Incident photon flux at Int(i)
   phi(i) = G*Int(i);
   
   
   %% Get the PL
   U = dfana.calcU(sol_OC(i));
        for j = 1:length(t)
            PLt(i,j) = trapz(x, U.btb(j,:));
        end
end


%% Get the ideality factor nid
for j = 1:length(t)
    nidt(:,j) = (par.q/(par.kB*par.T))*gradient(Voct(:,j),log(phi));
end

nidt_mean = mean(nidt,1);


figure(400)
hold on
for i = 1:length(sol_OC)
    plot(t, Voct)
end
hold off
xlabel('Time [s]')
ylabel('Voc [V]')

figure(401)
semilogx(t, nidt_mean)
xlabel('Time [s]')
ylabel('nid')
xlim([t(1), t(end)])
ylim([-1,4])

figure(402)
hold on
for i = 1:length(sol_OC)
    semilogy(t, PLt(i,:))
end
hold off
xlabel('Time [s]')
ylabel('PL [cm-2s-1]')

end