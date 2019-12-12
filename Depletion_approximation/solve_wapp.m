function [solwapp1, solwapp3] = solve_wapp
% Obtain expressions for wapp1 and Qapp3 in terms of Q

syms Q q N1 N2 N3 w1 w2 w3 eppx epp1 epp2 epp3 wapp1 wapp3 xf rho d0 d1 d2 d3 V wapp1 wapp3

% Potential at right hand boundary
eq1 = V == (q/(2*epp3))*(2*(N1.*(w1 + wapp1) - N3.*(w3 + wapp3)).*(d3 - d2) + N3.*(w3 + wapp3).^2) +...
                    (q/(2*epp2))*(2*(N1.*(w1 + wapp1)- N2*w2)*(d2-d1) +  2*N2.*w2^2) + q.*N1.*(w1 + wapp1).^2/(2*epp1);

% Field at right hand boundary
eq2 = 0 == q*(N3*(w3 + wapp3) - N1*(w1 + wapp1))/epp3;
                                        
eq2 = solve(eq2, wapp3);

eq3 = subs(eq1, wapp3, eq2);

solwapp1 = solve(eq3, wapp1);
solwapp3 = subs(eq3, wapp1, solwapp1);

solwapp1 = simplify(solwapp1, 'Steps', 100);
solwapp3 = simplify(solwapp3, 'Steps', 100);

end