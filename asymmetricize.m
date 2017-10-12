function [sol,newVoc] = asymmetricize(ssol, BC)
% break the symmetrical solution in two halves and stabilize the
% unsymmetrical model at pseudo-OC conditions applying Vapp equal to Voc
% taken from the symmetrical one

p = ssol.params;
p.OC = 0; % without setting OC to 0 BC gets ignored, BC 0 is needed for the asymmetric solution
p.BC = BC; % 

Voc = ssol.Voc(end); % get the Voc value in the middle of the symmetrical solution at the final time step
p.Vapp = Voc; % use potential value in the middle as new applied voltage, this is needed for BC 0
disp(Voc)

p.calcJ = 2; % verify that we're at open circuit equivalent conditions

centerIndex = ceil(length(ssol.x) / 2); % get the index of the middle of the x mesh
sol_ic.sol = ssol.sol(end, 1:centerIndex, :); % take the first spatial half of the matrix, at last time step
sol_ic.x = ssol.x(1:centerIndex); % providing just sol to pindrift is not enough, a more complete structure is needed, including at least x mesh

sol = pindrift(sol_ic, p); % should already be stable...
newVoc = sol.Efn(end) - sol.Efp(1);
disp(newVoc)

%assignin('base', 'sol', sol)
