function solve_interfaces

syms n x nl nr c1 c2 gam r alph mu kB T d p p0 bet

n = c1*exp(alph*x) - r*x/(alph*mu*kB*T) + c2;

% n = subs(n, c1, nl - (nr - gam*nl + r*d/(alph*mu*kB*T))/(1-gam));
% n = subs(n, c2, (nr - gam*nl + r*d/(alph*mu*kB*T))/(1-gam))
c1_express = (nl - nr - r*d/(alph*mu*kB*T))/(1-gam);
c2_express = nr - gam*c1 +r*d/(alph*mu*kB*T);
c2_express = subs(c2_express, c1, c1_express)

n = subs(n, c1, c1_express);
n = subs(n, c2, c2_express)

%n = subs(n, nr, gam*nl)
%n = exp(-alph*x)*(nl + (nr - gam*nl + (d*r)/(T*alph*kB*mu))/(gam - 1)) - (nr - gam*nl + (d*r)/(T*alph*kB*mu))/(gam - 1) - (r*x)/(T*alph*kB*mu)

%% Fluxes
dndx = diff(n)

jn = mu*kB*T*(alph*n - dndx);
jn = simplify(jn)

jn = -(exp(-alph*x)*(r*exp(alph*x) - 2*alph*d*r - gam*r*exp(alph*x) + alph*d*r*exp(alph*x) - alph*r*x*exp(alph*x) + alph*gam*r*x*exp(alph*x) + 2*T*alph^2*kB*mu*nl - 2*T*alph^2*kB*mu*nr + T*alph^2*kB*mu*nr*exp(alph*x) - T*alph^2*gam*kB*mu*nl*exp(alph*x)))/(alph*(gam - 1))

nr = gam*nl;
r = 0;

jn = subs(jn)

end

