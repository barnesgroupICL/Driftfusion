
%% Input parameters
n_maj = 1e18;
alpha = -1e8;
d = 1e-7;
points = 200;
mu = 10;

% Ideal case
r = 0;
n_min = n_maj*exp(alpha*d);
[x, n_ideal] = interface_analytical(n_maj, n_min, alpha, d, points, mu, r);

% Mobility limited
n_min = 1e15;
[x, n_mob_limited] = interface_analytical(n_maj, n_min, alpha, d, points, mu, r);

% Recombination limited
r = 1e30;
n_min = n_maj*exp(alpha*d);
[x, n_rec_limited] = interface_analytical(n_maj, n_min, alpha, d, points, mu, r);

%% Plot
figure(300)
semilogy(x, n_ideal, x, n_mob_limited, x, n_rec_limited)
ylabel('n (cm-3)')
xlabel('x (cm)')
ylim([1e12, 1e18])
legend('Ideal', 'Mob limited', 'Rec limited')