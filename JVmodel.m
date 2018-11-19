function JVmodel

Jph = -20e-3;
J0 = 2.9e11;
V = 0:0.01:1.4;

J = Jph + J0*exp(V/(2*0.0257)-1);

figure(600)
plot(V, J)
xlabel('V [V]')
ylabel('J [A]')
%ylim([-30e-3, 30e-3])
end