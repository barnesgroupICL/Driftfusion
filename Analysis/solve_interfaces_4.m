function solve_interfaces_4
% Solves for the interfacial carrier densities
% Boundary conditions: n(0) = ns, jn(0) = js
% Assumptions: constant rec rate, r

syms n xn xp ns ps jns jps c1 c2 gam r alph bet mue muh kB T d p p0 bet dndx ni taun taup nt pt

r_eqn = r == (ns*ps - ni^2)/(taun*(ps +pt) + taup*(ns + nt));

n_eqn = n == ns*exp(alph*xn) + jns/(alph*mue*kB*T)*(1-exp(alph*xn)) - (r/(alph^2 *mue*kB*T))*(1- exp(alph*xn) + alph*xn);

n_eqn = subs(n_eqn, r, rhs(r_eqn))
n_bc = subs(n_eqn, n, ns)
n_bc = subs(n_bc, xn, 0)

ps = solve(n_bc, ps)

p_eqn = p == ps*exp(bet*x) + jps/(bet*muh*kB*T)*(1-exp(bet*x)) - (r/(bet^2 *muh*kB*T))*(1- exp(bet*x) + bet*x);

p_eqn = subs(p_eqn, r, r_exp)
p_bc = subs(p_eqn, p, ps)
p_bc = subs(p_bc, x, 0)

end

