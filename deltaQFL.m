function deltaQFL = (n_integrated, p_integrated, Nc, Nv, d)
%% Universal constants
kB = 8.617330350e-5;     % Boltzmann constant [eV K^-1]
T = 300;
q = 1;

nbar = n_integrated/d;
pbar = p_integrated/d;

deltaQFL = (kB*T/q)*(ln(nbar/Nc) - ln(pbar.Nv));

end