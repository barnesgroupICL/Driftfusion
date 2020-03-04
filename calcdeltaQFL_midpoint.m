function [deltaQFL, np] = calcdeltaQFL_midpoint(sol, Nc, Nv, d, Eg)
%% Universal constants
kB = 8.617330350e-5;     % Boltzmann constant [eV K^-1]
T = 300;
q = 1;

[u,t,x,par,dev,n,p,a,c,V] = dfana.splitsol(sol);

p_midactive = find(x<=par.d_midactive);
p_midactive = p_midactive(end);

n = n(end, p_midactive);
p = p(end, p_midactive);
np = n*p;

deltaQFL = Eg + (kB*T/q)*(log(n/Nc) + log(p/Nv));

end