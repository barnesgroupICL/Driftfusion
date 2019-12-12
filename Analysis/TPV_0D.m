function TPV0Dstruct = TPV_0D(par, bias_int, pulse_int, tmax, duty)

%Universal constants
kBT = par.kB*par.T;        %Thermal energy eV
al = par.active_layer;

t_pulse = tmax*(duty*1e-2);
t_decay = tmax - t_pulse;

t1 = -t_pulse:1e-8:0;
t2 = 0:1e-8:t_decay;
t_tot = [t1, t2+t1(end)];

deltaG = pulse_int*par.g0(al);    % cm3s-1

k1 = 0;           % First order rate constant
k2 = par.krad(al);         % Second order rate constant

% Based on second order recombination only
n0 = ((bias_int*par.g0(al)/k2)+par.ni(al)^2)^0.5;

C = k1 + 2*k2*n0;      % Could go to higher orders- this constant falls out from discarding higher orders of deltan/

deltan = (deltaG/C)*(1-exp(-(t1+t_pulse)*C));

% TPV decay
deltan_fall = deltan(end)*exp(-t2*C);
deltan_tot = [deltan, deltan_fall];

% Intrinsic material
V0 = 2*kBT*log(n0/par.ni(al));
deltaV = 2*kBT*log((deltan_tot+n0)/par.ni(al)) - V0;
deltaVapprox = 2*kBT*deltan_tot/n0;

%deltaV2 = 2*kBT*log(1+(deltaG/(C*n0))*(1-exp(-C*t2)));

deltaV3 = 2*kBT*log(1+(deltan_tot/n0));

% Doped material
V0 = 2*kBT*log(n0/par.ni(al));
deltaV_dope = 2*kBT*log((deltan_tot+n0)/par.ni(al)) - V0;

% figure(101)
% plot(t_tot, deltan_tot)%, t, deltan_fall);
% xlabel('time [s]');
% ylabel('Charge density [cm-3]');

% Shift time
%t_tot = t_tot-t_pulse;

figure(600)
plot(t_tot, deltaV, 'k--');
xlabel('Time [s]');
ylabel('Change in Photovoltage, \Delta V [V]');

figure(601)
semilogy(t_tot, deltaV, 'k--');
xlabel('Time [s]');
ylabel('Change in Photovoltage, \Delta V [V]');

TPV0Dstruct.t = t_tot;
TPV0Dstruct.deltaVapprox = deltaVapprox;
TPV0Dstruct.V0 = V0;
TPV0Dstruct.g0 = par.g0(al);
TPV0Dstruct.n0 = n0;
TPV0Dstruct.deltan = deltan_tot;

end

