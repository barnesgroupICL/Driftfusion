function solQapp = solve_Qapp1
syms Q q N1 N2 N3 w1 w2 w3 eppx epp1 epp2 epp3 Qapp  xf rho d0 d1 d2 d3 V

% Potential solution at right hand boundary
eq1 = V == (q/(2*epp3))*(2*(N1.*(w1 + wapp1) - N3.*(w3 + wapp3)).*(x(ii)- d2) + N3.*(w3 + wapp3).^2) +...
                    (q/(2*epp2))*(2*(N1.*(w1 + wapp1)- N2*w2)*(d2-d1) +  2*N2.*w2^2) + q.*N1.*(w1 + wapp1).^2/(2*epp1);
            
                
                
eq1 = subs(eq1, wapp3, eq2);                
                
solQapp = solve(eq1, Qapp);

solQapp = simplify(solQapp, 'Steps', 100);

end