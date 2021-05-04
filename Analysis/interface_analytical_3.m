function [x, n, jn] = interface_analytical_3(ns, js, alpha, d, points, mu, r)

T = 300;
kB = 0.0257/300;
kT = kB*T;
alph = alpha;

x = 0:d/points:d;

n = ns*exp(alpha*x) +(js/(alpha*mu*kB*T))*(1-exp(alpha*x)) - (r/(alpha^2*mu*kB*T))*(1-exp(alpha*x) + alpha*x);
jn = js - r*x;

figure(3002)
plot(x, 1-exp(alpha*x), x, alpha*x)
ylabel('components of r')
xlabel('x')
legend('1-exp(alpha*x)', 'alpha*x')
% rec_term1 = exp(alph*x).*(r./alph)./(T*alph*kB*mu);
% rec_term2 = ones(1,length(x)) * - (r./(T*alph^2*kB*mu));
% rec_term3 = - (r*x)./(T*alph*kB*mu);
% 
% figure(301)
% plot(x, rec_term1, x, rec_term2, x, rec_term3)
% xlabel('x')
% ylabel('r')
% legend('r1', 'r2', 'r3')
end