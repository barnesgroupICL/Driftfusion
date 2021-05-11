function solve_interfaces_6
% Solves for the interfacial carrier densities
% Boundary conditions: n(0) = ns, jn(0) = js
% Assumptions: constant rec rate, r
% This r is equated in the two sets of eqns

syms q kso n x xn xp ns ps jn jp jns jps c1 c2 gam r rn rp alph bet mue muh kB T d p p0 bet dndx ni taun taup nt pt

%% SRH
r_eqn = r == (ns*ps - ni^2)/(taun*(ps +pt) + taup*(ns + nt));

%% Constant mobility
n_eqn = n == ns*exp(alph*xn) + (q*jns/(alph*mue*kB*T))*(1-exp(alph*xn)) - (rn/(alph^2 *mue*kB*T))*(1- exp(alph*xn) + alph*xn);
p_eqn = p == ps*exp(bet*xp) + (q*jps/(bet*muh*kB*T))*(1-exp(bet*xp)) - (rp/(bet^2 *mue*kB*T))*(1- exp(bet*xp) + bet*xp);

subs(n_eqn, rn, (jns-jn)/xn)
simplify(n_eqn)
subs(p_eqn, rp, (jps-jp)/xn)
simplify(p_eqn)


%% exponential mobility
%n_eqn = n == (ns - q*jns*x/(mue*kB*T) + (r*x^2)/(2*kB*T))*exp(alph*x)

n_eqn = subs(n_eqn, r, rhs(r_eqn))
ns_expr = solve(n_eqn, ns)
ns_expr = simplify(ns_expr)

%% Constant mobility

%% exponential mobility
%p_eqn = p == (ps - q*jps*x/(muh*kB*T) + (r*x^2)/(2*kB*T))*exp(bet*x)

p_eqn = subs(p_eqn, r, rhs(r_eqn))
ps_expr = solve(p_eqn, ps)
ps_expr = simplify(ps_expr)

ps_expr = subs(ps_expr, ns, ns_expr);
ps_expr = subs(ps_expr, xn, d-xp);
ps_expr = simplify(ps_expr,100)
end

