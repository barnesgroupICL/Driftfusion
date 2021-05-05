function solve_interfaces_4
% Solves for the interfacial carrier densities
% Boundary conditions: n(0) = ns, jn(0) = js
% Assumptions: constant rec rate, r

syms q kso n x ns ps jns jps c1 c2 gam r alph bet mue muh kB T d p p0 bet dndx ni taun taup nt pt

%% SRH
r_eqn = r == (ns*ps - ni^2)/(taun*(ps +pt) + taup*(ns + nt));

%% Second order
%r_eqn = r == kso*ns*ps;

%% Constant mobility
n_eqn = n == ns*exp(alph*x) + q*jns/(alph*mue*kB*T)*(1-exp(alph*x)) - (r/(alph^2 *mue*kB*T))*(1- exp(alph*x) + alph*x);

%% exponential mobility
%n_eqn = n == (ns - q*jns*x/(mue*kB*T) + (r*x^2)/(2*kB*T))*exp(alph*x)

n_eqn = subs(n_eqn, r, rhs(r_eqn))
ns_expr = solve(n_eqn, ns)
ns_expr = simplify(ns_expr)

%% Constant mobility
p_eqn = p == ps*exp(alph*x) + q*jps/(alph*mue*kB*T)*(1-exp(alph*x)) - (r/(alph^2 *mue*kB*T))*(1- exp(alph*x) + alph*x);

%% exponential mobility
%p_eqn = p == (ps - q*jps*x/(muh*kB*T) + (r*x^2)/(2*kB*T))*exp(bet*x)

p_eqn = subs(p_eqn, r, rhs(r_eqn))
ps_expr = solve(p_eqn, ps)
ps_expr = simplify(ps_expr)

ps_expr = subs(ps_expr, ns, ns_expr);
ps_expr = simplify(ps_expr)
end

