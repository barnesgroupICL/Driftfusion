function SCLCAna(solstruct)

% Simple structure names
sol = solstruct.sol;
p = solstruct.p;
Vapp_arr = solstruct.Vapp;

%% ANALYSIS %%
% split the solution into its component parts (e.g. electrons, holes and efield)
n = sol(:,:,1);     % electrons
P = sol(:,:,2);     % holes
a = sol(:,:,3);     % mobile ions
V = sol(:,:,4);     % electric potential

%% Currents from right-hand boundary

Jn_r = p.sn_r*(n(:, end) - p.n0n)*-p.e;
Jp_r = p.sp_r*(P(:, end) - p.p0n)*p.e;
Jtot = Jn_r + Jp_r;

gradJV = gradient(log(Jtot), log(Vapp_arr));

% Mobility calculation based on Mott-Gurney law
mu_MG = gradient(Jtot, (Vapp_arr.^2))*(8/9)*(p.ti^3)/(p.eppi*p.e);

figure(200)
loglog(Vapp_arr, Jtot)
xlabel('Vapp [V]')
ylabel('Current Density [A cm^-2]');
grid off;

figure(201)
semilogx(Vapp_arr, gradJV)
xlabel('V_{app} [V]')
ylabel('dJ/dV')

figure(202)
semilogx(Vapp_arr, mu_MG)
xlabel('Vapp [V]')
ylabel('Calculated mobility')

end