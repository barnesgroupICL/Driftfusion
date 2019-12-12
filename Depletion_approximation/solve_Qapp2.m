function solQapp = solve_Qapp2
syms Q q N1 N2 N3 w1 w2 w3 eppx epp1 epp2 epp3 Qapp rho d0 d1 d2 d3 V wapp1 wapp3

% Potential solution at right hand boundary
eq1 = V == (q/(2*epp3))*(N3.*(w3 + wapp3).^2) + (q/(2*epp2))*(2*(N1.*(w1 + wapp1)- N2*w2)*(d2-d1) +...
                2*N2.*w2^2) + q.*N1.*(w1 + wapp1).^2/(2*epp1);

exp1 = w1 == Q/(q*N1);
exp2 = w2 == Q/(q*N2);       
exp3 = w3 == Q/(q*N3);  
 
exp4 = wapp1 == Qapp/(q*N1);
exp5 = wapp3 == Qapp/(q*N3);
                
eq1 = subs(eq1, w1, rhs(exp1));                
eq1 = subs(eq1, w2, rhs(exp2));
eq1 = subs(eq1, w3, rhs(exp3));
eq1 = subs(eq1, wapp1, rhs(exp4));                
eq1 = subs(eq1, wapp3, rhs(exp5));

eq1_simple = V == (N1*q*(Q/(N1*q) + Qapp/(N1*q))^2)/(2*epp1) - (q*((2*N1*(Q/(N1*q) + Qapp/(N1*q)) - (2*Q)/q)*(d1 - d2) - (2*Q^2)/(N2*q^2)))/(2*epp2) + (N3*q*(Q/(N3*q) + Qapp/(N3*q))^2)/(2*epp3);

solQapp = solve(eq1, Qapp);

solQapp = simplify(solQapp, 'Steps', 100);

end