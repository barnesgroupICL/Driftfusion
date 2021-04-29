function [x, n] = interface_analytical(n_maj, n_min, alpha, d, points, mu, r)

nr = n_min;
nl = n_maj;
alph = alpha;
gam = exp(alpha*d);

T = 300;
kB = 0.0257/300;
kT = kB*T;

x = 0:d/points:d;

%n = n_maj*exp(alpha*x) + ((n_min - n_maj*exp(alpha*d) + r*d/(mu*kT*alpha))*(exp(alpha*x) - 1)/(exp(alpha*d) - 1)) - r*x/(mu*kT*alpha);

n = nr + (exp(alph*x)*(nr - nl + (d*r)/(T*alph*kB*mu)))/(gam - 1) - (gam*(nr - nl + (d*r)/(T*alph*kB*mu)))/(gam - 1) + (d*r)/(T*alph*kB*mu) - (r*x)/(T*alph*kB*mu)

end