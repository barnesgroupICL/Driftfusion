function transistormodel(Vn, Js1, Js2)

% For the mixed regime
e = 1.602e-19;
kT = 0.0257;
Jph = -20e-3;
V = 0:0.01:1.4;
Vn = Vn*ones(1, length(V));
Vrec = V-Vn;
Vcol = -Vn; 

J1 = Js1*(exp(V/kT)-1);
J2 = Js2*(exp(V/kT)-1);

J = Jph + J1;

figure(600)
semilogy(V, J1, V, J2)
xlabel('V [V]')
ylabel('J [A]')
legend('J1', 'J2')
%ylim([-30e-3, 30e-3])

figure(601)
plot(V, J)
xlabel('V [V]')
ylabel('J [A]')
%legend('J1', 'J2')
ylim([-30e-3, 30e-3])

end