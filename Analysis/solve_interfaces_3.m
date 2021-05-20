function solve_interfaces_3
% Solves for the interfacial carrier densities
% Boundary conditions 1: n(0) = ns, jn(0) = js
% Boundary conditions 2: n(d) = ns, jn(0) = js
% Assumptions: 
% - constant rec rate, r
% - constant mobility
% - constant alpha
% - constant field

syms n dndx d2ndx2 x d ns jn js c1 c2 gam r alph mu mus kB T p p0 bet 

neqn = n == (c1/alph)*exp(alph*x) - r*x/(alph*mu*kB*T) + c2;
dndxeqn = dndx == diff(rhs(neqn), x);
d2ndx2eqn = d2ndx2 == diff(rhs(dndxeqn), x)

jn_eqn = jn == mu*kB*T*(alph*n - dndx)
jn_eqn = subs(jn_eqn, n, rhs(neqn))
jn_eqn = subs(jn_eqn, dndx, rhs(dndxeqn))

%% Boundary conditions
%% n(0) = ns, jn(0) = js
neqn_bc1 = neqn;
neqn_bc1 = subs(neqn, x, 0);
neqn_bc1 = subs(neqn_bc1, n, ns);
c1sol = solve(neqn_bc1, c1)

jn_eqn_bc1 = jn_eqn;
jn_eqn_bc1 = subs(jn_eqn, x, 0);
jn_eqn_bc1 = subs(jn_eqn_bc1, jn, js);
jn_eqn_bc1 = simplify(jn_eqn_bc1);
c2sol = solve(jn_eqn_bc1, c2);

neqn_bc1_final = subs(neqn, c1, c1sol)
neqn_bc1_final = subs(neqn_bc1_final, c2, c2sol)

d2ndx2eqn = d2ndx2 == diff(rhs(neqn_bc1_final), x)

%% n(d) = ns, jn(0) = js
neqn_bc2 = neqn;
neqn_bc2 = subs(neqn, x, d)
neqn_bc2 = subs(neqn_bc2, n, ns)
c1sol = solve(neqn_bc2, c1)

jn_eqn_bc2 = jn_eqn;
jn_eqn_bc2 = subs(jn_eqn, x, 0)
jn_eqn_bc2 = subs(jn_eqn_bc2, jn, js)
jn_eqn_bc2 = subs(jn_eqn_bc2, c1, c1sol)
c2sol = solve(jn_eqn_bc2, c2)

neqn_bc2_final = subs(neqn, c1, c1sol)
neqn_bc2_final = subs(neqn_bc2_final, c2, c2sol)
neqn_bc2_final = subs(neqn_bc2_final, r, (js-jn)/x)

end

