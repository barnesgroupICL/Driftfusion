function asymstruct = asymmetricize(symstruct)
%ASYMMETRICIZE - Break a symmetrical solution in two halves
% Please note that further stabilization could be required
%
% Syntax:  asymstruct = asymmetricize(symstruct)
%
% Inputs:
%   SYMSTRUCT - a symmetric struct as created by PINDRIFT using the OC
%     parameter set to 1
%
% Outputs:
%   ASYMSTRUCT - an asymmetric struct as created by PINDRIFT using the OC
%     parameter set to 0 and applied voltage identical to the open circuit
%     voltage taken from the input SYMSTRUCT
%
% Example:
%   sol_i_light_OC = asymmetricize(ssol_i_SR_1S)
%     take the first half of symmetrical solution at open circuit
%
% Other m-files required: pindrift
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift.

% Author: Ilario Gelmetti, Ph.D. student, perovskite photovoltaics
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Phil Calado, Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% October 2017; Last revision: May 2018

%------------- BEGIN CODE --------------

p = symstruct.p;
p.OC = 0; % OC 0 is needed for the asymmetric solution
p.tpoints = 10; % rough, just for re-stabilization at Vapp
p.Ana = 0;

% get the Voc value at the final time step
[Voc, ~, ~] = pinana(symstruct);
p.Vapp = Voc(end); % use potential value in the middle as new applied voltage

% set an initial time for stabilization tmax
if p.muion
    p.tmax = min(1e0, 2^(-log10(p.muion)) / 10 + 2^(-log10(p.mue(1))));
else
    p.tmax = min(1e-3, 2^(-log10(p.mue(1))));
end

%% halve the solution
centerIndex = ceil(length(symstruct.x) / 2); % get the index of the middle of the x mesh
sol_ic.sol = symstruct.sol(end, 1:centerIndex, :); % take the first spatial half of the matrix, at last time step
sol_ic.x = symstruct.x(1:centerIndex); % providing just sol to pindrift is not enough, a more complete structure is needed, including at least x mesh

%% stabilize the divided solution

p.BC = 1; % safe boundary conditions for breaking solutions
asymstruct = pindrift(sol_ic, p); % should already be stable, this run of pindrift should be fast...
if symstruct.p.BC ~= 1 % in case the boundary conditions set was not 1, go to the correct BC
    p.BC = symstruct.p.BC; % re-establish the original BC
    asymstruct = pindrift(asymstruct, p);
end

%------------- END OF CODE --------------
