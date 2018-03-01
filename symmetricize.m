function symsol = symmetricize(sol)
% Function to symmetricize solution

x = sol.x;

% Symmetricise main solution matrix
sol1 = sol.sol(:, :, :);
sol2 = fliplr(sol1(:, :, :));
sol2 = sol2(:, 2:end, :);               % Delete middle point
symsol.sol = [sol1, sol2];              % Concatenate

% Flip xmesh
x1 = x;
x2 = 2*x(end) - fliplr(x);
x2 = x2(2:end);    
symsol.x = [x1, x2];

end
