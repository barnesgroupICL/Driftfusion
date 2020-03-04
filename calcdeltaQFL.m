function [deltaQFL, nbarpbar] = calcdeltaQFL(n_integrated, p_integrated, Nc, Nv, d, Eg)
%% Universal constants
kB = 8.617330350e-5;     % Boltzmann constant [eV K^-1]
T = 300;
q = 1;

nbar = n_integrated/d;
pbar = p_integrated/d;
nbarpbar = nbar*pbar;

deltaQFL = Eg + (kB*T/q)*(log(nbar/Nc) + log(pbar/Nv));

end