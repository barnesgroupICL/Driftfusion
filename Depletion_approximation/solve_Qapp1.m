function solQapp = solve_Qapp1
syms Q q N1 N2 N3 w1 w2 w3 eppx epp1 epp2 epp3 Qapp  xf rho d0 d1 d2 d3 V

% Potential solution at right hand boundary
eq1 = V == (1/2)*(((Q+Qapp)^2)*(1/(q*N1*epp1) + 1/(q*N3*epp3))) + Qapp*(d2-d1)/epp2 + Q^2/(q*N2*epp2);
               
solQapp = solve(eq1, Qapp);

solQapp = simplify(solQapp, 'Steps', 100);

end