function solve_interfaces_4
% Solves for the interfacial carrier densities
% Boundary conditions: n(0) = ns, jn(0) = js
% Assumptions: constant rec rate, r

syms n x ns jn js c1 c2 gam r alph mu mus kB T d p p0 bet 

neq = (n == (c1/alph)*exp(alph*x) - r*x/(alph*mu*kB*T) + c2);
% mu_expr = mus*exp(-alph*x);
% neq = subs(neq, mu, mu_expr)

dndx = diff(rhs(neq), x)

jneq = (jn == mu*kB*T*(alph*rhs(neq) - dndx))
jn = simplify(rhs(jneq))
djndx = diff(rhs(jneq))

c1_expr = solve(neq, c1)
c2_expr = solve(jn, c2)

neq = subs(neq, c1, c1_expr)
neq = subs(neq, c2, c2_expr)

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

