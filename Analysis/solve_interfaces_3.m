function solve_interfaces_3

syms n x ns js c1 c2 gam r alph mu kB T d p p0 bet


jn = mu*kB*T*(alph*n - diff(n));
djndx = diff(jn)

solve(diff(jn) + r = 0, n)


% c1_express = ns - c2
% c2_express = (1/(alph*mu*kB*T))*(js - r*(1/alph))
% c1_express = subs(c1_express, c2, c2_express)
% 
% n = subs(n, c1, c1_express);
% n = subs(n, c2, c2_express)
% 
% 
% %% Fluxes
% dndx = diff(n)
% 
% jn = mu*kB*T*(alph*n - dndx);
% jn = simplify(jn)
%n = subs(n, nr, gam*nl)
%n = exp(-alph*x)*(nl + (nr - gam*nl + (d*r)/(T*alph*kB*mu))/(gam - 1)) - (nr - gam*nl + (d*r)/(T*alph*kB*mu))/(gam - 1) - (r*x)/(T*alph*kB*mu)


end

