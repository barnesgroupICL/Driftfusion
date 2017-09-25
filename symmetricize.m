function symsol = symmetricize(sol)

% Function to symmetricize solution

xp = length(sol.x);

x = sol.x
n = sol.n
p = sol.p
V = sol.V

% First copy solution
symsol = sol;

% Symmetricise main solution matrix
sol1(:, 1:xp, :) = symsol.sol(:,:,:);
sol2 = fliplr(sol1(:,:,:));
sol2 = sol2(:, 2:end, :);               % Delete middle point
symsol.sol = [sol1, sol2];              % Concatenate

% Flip xmesh
x1 = x;
x2 = 2*x(end) - fliplr(x);
x2 = x2(2:end);    
symsol.x = [x1, x2]

% symsol.n(end, 1:xp) = sol.n;
% symsol.n(end, xp+1:2*xp) = fliplr(sol.n(end, :));
% 
% symsol.p(end, 1:xp) = sol.p;
% symsol.p(end, xp+1:2*xp) = fliplr(sol.p(end, :));
% 
% symsol.i(end, 1:xp) = sol.i;
% symsol.i(end, xp+1:2*xp) = fliplr(sol.i(end, :));
% 
% symsol.V(end, 1:xp) = sol.V;
% symsol.V(end, xp+1:2*xp) = fliplr(sol.V(end, :));

assignin('base', 'symsol', symsol);

end