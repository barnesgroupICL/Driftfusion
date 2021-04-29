
%% Input parameters
n_maj = 1e18;
alpha = -1e8;
d = 1e-7;
points = 400;

%% Ideal case
js = -1e16;
mu = 1e6;
r = 0;
[x, n_ideal, jn_ideal] = interface_analytical_2(n_maj, js, alpha, d, points, mu, r);

%% Mobility limited
js = -1e16;
mu = 1e-6;
r = 0;
[x, n_mob_limited, jn_mob_limited] = interface_analytical_2(n_maj, js, alpha, d, points, mu, r);

%% Recombination limited
js = 1e23;  % This has to change sign as carriers are no longer being extracted, but sucked into the recombination zone
mu = 10;
r = 1e30;
[x, n_rec_limited, jn_rec_limited] = interface_analytical_2(n_maj, js, alpha, d, points, mu, r);

%% Plot
figure(300)
semilogy(x, n_ideal, x, n_mob_limited, x, n_rec_limited)
ylabel('n (cm-3)')
xlabel('x (cm)')
%ylim([1e12, 1e18])
legend('Ideal', 'Mob limited', 'Rec limited')

figure(3001)
plot(x, jn_ideal, x, jn_mob_limited, x, jn_rec_limited)
ylabel('jn (cm-2s-1)')
xlabel('x (cm)')
%ylim([1e12, 1e18])
legend('Ideal', 'Mob limited', 'Rec limited')