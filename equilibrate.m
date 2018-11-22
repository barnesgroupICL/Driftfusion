function [sol_eq, sol_eq_SR, sol_i_eq, sol_i_eq_SR, ssol_eq, ssol_eq_SR, ssol_i_eq, ssol_i_eq_SR, sol_1S, sol_1S_SR, sol_i_1S, sol_i_1S_SR, ssol_1S, ssol_1S_SR, ssol_i_1S, ssol_i_1S_SR] = equilibrate(params)
%EQUILIBRATE Uses analytical initial conditions and runs to equilibrium and steady state
% Takes the parameters from pinParams.m file and tries
% to obtain an equilibrium solution (as if the device has been left for
% a long period of time). This solution can then be used as accurate
% initial conditions for other simulations, e.g. a JV scan.
% Note that tmax is consistently adjusted to appropriate values for to
% ensure there are numerous mesh points where large gradients in the time
% dimension are present.
%
% Syntax:  [sol_eq, sol_eq_SR, sol_i_eq, sol_i_eq_SR, ssol_eq, ssol_eq_SR, ssol_i_eq, ssol_i_eq_SR, sol_light, sol_light_SR, sol_i_light, sol_i_light_SR, ssol_light, ssol_light_SR, ssol_i_light, ssol_i_light_SR] = EQUILIBRATE(PARAMS)
%
% Inputs:
%   PARAMS - optional, struct containing the needed parameters as obtained
%     from pinParams.m
%
% Outputs:
%   sol_eq - short circuit, dark, no mobile ionic defects, no SRH
%   sol_eq_SR - short circuit, dark, no mobile ionic defects, with SRH
%   sol_i_eq - short circuit, dark, mobile ionic defects, no SRH
%   sol_i_eq_SR - short circuit, dark, mobile ionic defects, with SRH
%   ssol_eq - open circuit, dark, no mobile ionic defects, no SRH
%   ssol_eq_SR - open circuit, dark, no mobile ionic defects, with SRH
%   ssol_i_eq - open circuit, dark, mobile ionic defects, no SRH
%   ssol_i_eq_SR - open circuit, dark, mobile ionic defects, with SRH
%   sol_1S - short circuit, 1 sun, no mobile ionic defects, no SRH
%   sol_1S_SR - short circuit, 1 sun, no mobile ionic defects, with SRH
%   sol_i_1S - short circuit, 1 sun, mobile ionic defects, no SRH
%   sol_i_1S_SR - short circuit, 1 sun, mobile ionic defects, with SRH
%   ssol_1S - open circuit, 1 sun, no mobile ionic defects, no SRH
%   ssol_1S_SR - open circuit, 1 sun, no mobile ionic defects, with SRH
%   ssol_i_1S - open circuit, 1 sun, mobile ionic defects, no SRH
%   ssol_i_1S_SR - open circuit, 1 sun, mobile ionic defects, with SRH
%
% Example:
%   [sol_eq, sol_eq_SR, sol_i_eq, sol_i_eq_SR, ssol_eq, ssol_eq_SR, ssol_i_eq,...
%       ssol_i_eq_SR, sol_1S, sol_1S_SR, sol_i_1S, sol_i_1S_SR,...
%       ssol_1S, ssol_1S_SR, ssol_i_1S, ssol_i_1S_SR] = equilibrate()
%     generate stabilized solutions using default parameters from
%     pinParams.m
%
% Other m-files required: pindrift, pinParams, stabilize
% Subfunctions: none
% MAT-files required: none
%
% See also pindrift, paramsStruct, equilibrate_minimal.

% Author: Phil Calado, Ph.D.
% Imperial College London
% Research Group Prof. Jenny Nelson
% email address: p.calado13@imperial.ac.uk
% Contributors: Ilario Gelmetti, Ph.D. student
% Institute of Chemical Research of Catalonia (ICIQ)
% Research Group Prof. Emilio Palomares
% email address: iochesonome@gmail.com
% Supervised by: Dr. Piers Barnes, Prof. Jenny Nelson
% Imperial College London
% 2015; Last revision: May 2018

%------------- BEGIN CODE --------------

tic;    % Start stopwatch

%% Initial arguments
% Setting sol.sol = 0 enables a parameters structure to be read into
% pindrift but indicates that the initial conditions should be the
% analytical solutions
sol0.sol = 0;

% if a params struct has been provided in input, use it instead of the
% pinParams file
if nargin
    p = params;
else
    p = pinParams;
end

% Store initial parameters
original_p = p;

%% Start with low recombination coefficients
p.klin = 0;
p.klincon = 0;
p.taun_etl = 1e6; % [s] SRH time constant for electrons in ETL
p.taup_etl = 1e6; % [s] SRH time constant for holes in ETL
p.taun_htl = 1e6; % [s] SRH time constant for electrons in HTL
p.taup_htl = 1e6; % [s] SRH time constant for holes in HTL
p.taun_i = 1e6;
p.taup_i = 1e6;

% Raditative recombination could also be set to low values initially if required. 
% p.krad = 1e-20;
% p.kradetl = 1e-20;
% p.kradhtl = 1e-20;

%% General initial parameters

% likely we're just interested in the last time point, so simulating just
% a few is enough
p.tpoints = 20;

p.Ana = 0; % analysis is not needed
p.JV = 0;
p.Vapp = 0;
p.Int = 0;
p.pulseon = 0; 
p.OC = 0;
p.tmesh_type = 2;
p.tmax = 1e-9;
p.t0 = p.tmax / 1e4;
p.figson = 0; % reduce annoyance of figures popping up
p.calcJ = 0; % current calculation is not needed

%% Switch off mobilities
p.mue_i = 0; % electron mobility in intrinsic
p.muh_i = 0; % hole mobility in intrinsic
p.mue_p = 0; % electron mobility in p-type
p.muh_p = 0; % hole mobility in n-type
p.mue_n = 0; % electron mobility in p-type
p.muh_n = 0; % hole mobility in n-type
p.mui = 0; % ionic mobility in intrinsic

% Switch off extraction and recombination, used in BC 3
p.sn_ext = 0;
p.sn_rec = 0;
p.sp_ext = 0;
p.sp_rec = 0;

%% Initial solution with zero mobility
disp('Initial solution, zero mobility')
sol1 = pindrift(sol0, p);

%% Mobility with mobility switched on

% switch on electron and hole mobility
p.mue_i = original_p.mue_i; % electron mobility in intrinsic
p.muh_i = original_p.muh_i; % hole mobility in intrinsic
p.mue_p = original_p.mue_p; % electron mobility in p-type
p.muh_p = original_p.muh_p; % hole mobility in n-type
p.mue_n = original_p.mue_n; % electron mobility in p-type
p.muh_n = original_p.muh_n; % hole mobility in n-type

% Switch on extraction and recombination, used in BC 3
p.sn_ext = original_p.sn_ext;
p.sn_rec = original_p.sn_rec;
p.sp_ext = original_p.sp_ext;
p.sp_rec = original_p.sp_rec;

disp('Solution with mobility switched on')
sol2 = pindrift(sol1, p);

p.tmax = 1e-3;
p.t0 = p.tmax / 1e8;

sol_eq = pindrift(sol2, p);
sol_eq = stabilize(sol_eq); % go to steady state

sol_eq_p = p; % temporarily save params

%% Set up solution for open circuit
disp('Switching boundary conditions to zero flux')
p.BC = 0;
p.tmax = 1e-9;
p.t0 = p.tmax / 1e3;

sol3 = pindrift(sol_eq, p);

%% Symmetricise the solution
disp('Symmetricise solution for open circuit')
symsol = symmetricize(sol3);

%% Equilibrium solution with mirrored cell and OC boundary conditions, mobility zero
disp('Initial equilibrium open circuit solution')

% update p to symmetric solution (x, xpoints, OC)
p = symsol.p;

p.BC = original_p.BC;

% Switch off mobilities
p.mue_i = 0; % electron mobility in intrinsic
p.muh_i = 0; % hole mobility in intrinsic
p.mue_p = 0; % electron mobility in p-type
p.muh_p = 0; % hole mobility in n-type
p.mue_n = 0; % electron mobility in p-type
p.muh_n = 0; % hole mobility in n-type
p.mui = 0; % ionic mobility in intrinsic

ssol = pindrift(symsol, p);

%% OC with mobility switched on
disp('Open circuit solution with mobility switched on')
p.tmax = 1e-6;
p.t0 = p.tmax/1e3;

% switch on electron and hole mobility
p.mue_i = original_p.mue_i; % electron mobility in intrinsic
p.muh_i = original_p.muh_i; % hole mobility in intrinsic
p.mue_p = original_p.mue_p; % electron mobility in p-type
p.muh_p = original_p.muh_p; % hole mobility in n-type
p.mue_n = original_p.mue_n; % electron mobility in p-type
p.muh_n = original_p.muh_n; % hole mobility in n-type

ssol_eq = pindrift(ssol, p);
ssol_eq = stabilize(ssol_eq); % go to steady state

ssol_eq_p = p; % temporarily save params

%% Equilibrium solutions with ion mobility switched on - Closed circuit conditions
disp('Closed circuit equilibrium with ions')

p.OC = 0;
p.tmax = 1e-6;
p.t0 = p.tmax / 1e3;

% Ions are set to the correct speed indicated in pinParams
p.mui = original_p.mui;

sol_i_eq = pindrift(sol_eq, p);
sol_i_eq = stabilize(sol_i_eq); % go to steady state

sol_i_eq_p = p; % temporarily save params

%% Ion equilibrium with surface recombination
disp('Switching on surface recombination')
p.taun_etl = original_p.taun_etl;
p.taup_etl = original_p.taup_etl;
p.taun_htl = original_p.taun_htl;
p.taup_htl = original_p.taup_htl;
p.taun_i = original_p.taun_i;
p.taup_i = original_p.taup_i;

p.tmax = 1e-6;
p.t0 = p.tmax / 1e3;

sol_i_eq_SR = pindrift(sol_i_eq, p);
sol_i_eq_SR = stabilize(sol_i_eq_SR);

sol_i_eq_SR_p = p; % temporarily save params

% Switch off SR
p.taun_etl = 1e6;
p.taup_etl = 1e6;
p.taun_htl = 1e6;
p.taup_htl = 1e6;
p.taun_i = 1e6;
p.taup_i = 1e6;

%% Symmetricise closed circuit condition
disp('Symmetricise equilibriumion solution')
symsol = symmetricize(sol_i_eq);

% update p to symmetric solution (x, xpoints, OC)
p = symsol.p;

p.tmax = 1e-9;
p.t0 = p.tmax / 1e3;

% Switch off mobilities
p.mue_i = 0; % electron mobility in intrinsic
p.muh_i = 0; % hole mobility in intrinsic
p.mue_p = 0; % electron mobility in p-type
p.muh_p = 0; % hole mobility in n-type
p.mue_n = 0; % electron mobility in p-type
p.muh_n = 0; % hole mobility in n-type
p.mui = 0; % ionic mobility in intrinsic

%% OC condition with ions at equilbirium
disp('Open circuit equilibrium with ions')
ssol = pindrift(symsol, p);

% switch on electron and hole mobility
p.mue_i = original_p.mue_i; % electron mobility in intrinsic
p.muh_i = original_p.muh_i; % hole mobility in intrinsic
p.mue_p = original_p.mue_p; % electron mobility in p-type
p.muh_p = original_p.muh_p; % hole mobility in n-type
p.mue_n = original_p.mue_n; % electron mobility in p-type
p.muh_n = original_p.muh_n; % hole mobility in n-type

ssol = pindrift(ssol, p);

% Switch on ion mobility
p.tmax = 1e-6;
p.t0 = p.tmax / 1e3;
p.mui = original_p.mui; % this requires mui to be set in the original params

ssol_i_eq = pindrift(ssol, p);
ssol_i_eq = stabilize(ssol_i_eq);

ssol_i_eq_p = p; % temporarily save params

%% Dark, short circuit, surface recombination
disp("Dark, short circuit, surface recombination")
p = sol_eq_p;
p.taun_etl = original_p.taun_etl;
p.taup_etl = original_p.taup_etl;
p.taun_htl = original_p.taun_htl;
p.taup_htl = original_p.taup_htl;
p.taun_i = original_p.taun_i;
p.taup_i = original_p.taup_i;

sol_eq_SR = pindrift(sol_eq, p);
sol_eq_SR = stabilize(sol_eq_SR);

sol_eq_SR_p = p; % temporarily save params

%% Dark, open circuit, surface recombination
disp("Dark, open circuit, surface recombination")
p = ssol_eq_p;
p.taun_etl = original_p.taun_etl;
p.taup_etl = original_p.taup_etl;
p.taun_htl = original_p.taun_htl;
p.taup_htl = original_p.taup_htl;
p.taun_i = original_p.taun_i;
p.taup_i = original_p.taup_i;

ssol_eq_SR = pindrift(ssol_eq, p);
ssol_eq_SR = stabilize(ssol_eq_SR);

ssol_eq_SR_p = p; % temporarily save params

%% Dark, mobile ions, open circuit, surface recombination
disp("Dark, mobile ions, open circuit, surface recombination")
p = ssol_i_eq_p;
p.taun_etl = original_p.taun_etl;
p.taup_etl = original_p.taup_etl;
p.taun_htl = original_p.taun_htl;
p.taup_htl = original_p.taup_htl;
p.taun_i = original_p.taun_i;
p.taup_i = original_p.taup_i;

ssol_i_eq_SR = pindrift(ssol_i_eq, p);
ssol_i_eq_SR = stabilize(ssol_i_eq_SR);

ssol_i_eq_SR_p = p; % temporarily save params

%% Illuminated, short circuit
disp("Illuminated, short circuit")
p = sol_eq_p;
p.Int = original_p.Int;

sol_1S = pindrift(sol_eq, p);
sol_1S = stabilize(sol_1S);

%% Illuminated, short circuit, surface recombination
disp("Illuminated, short circuit, surface recombination")
p = sol_eq_SR_p;
p.Int = original_p.Int;

sol_1S_SR = pindrift(sol_eq_SR, p);
sol_1S_SR = stabilize(sol_1S_SR);

%% Illuminated, mobile ions, short circuit
disp("Illuminated, mobile ions, short circuit")
p = sol_i_eq_p;
p.Int = original_p.Int;

sol_i_1S = pindrift(sol_i_eq, p);
sol_i_1S = stabilize(sol_i_1S);

%% Illuminated, mobile ions, short circuit, surface recombination
disp("Illuminated, mobile ions, short circuit, surface recombination")
p = sol_i_eq_SR_p;
p.Int = original_p.Int;

sol_i_1S_SR = pindrift(sol_i_eq_SR, p);
sol_i_1S_SR = stabilize(sol_i_1S_SR);

%% Illuminated, open circuit
disp("Illuminated, open circuit")
p = ssol_eq_p;
p.Int = original_p.Int;
p.tmax = 1e-1;

ssol_1S = pindrift(ssol_eq, p);
ssol_1S = stabilize(ssol_1S);

%% Illuminated, open circuit, surface recombination
disp("Illuminated, open circuit, surface recombination")
p = ssol_eq_SR_p;
p.Int = original_p.Int;

ssol_1S_SR = pindrift(ssol_eq_SR, p);
ssol_1S_SR = stabilize(ssol_1S_SR);

%% Illuminated, mobile ions, open circuit
disp("Illuminated, mobile ions, open circuit")
p = ssol_i_eq_p;
p.Int = original_p.Int;

ssol_i_1S = pindrift(ssol_i_eq, p);
ssol_i_1S = stabilize(ssol_i_1S);

%% Illuminated, mobile ions, open circuit, surface recombination
disp("Illuminated, mobile ions, open circuit, surface recombination")
p = ssol_i_eq_SR_p;
p.Int = original_p.Int;

p.tmax = p.tmax / 1e1;
ssol_i_1S_SR = pindrift(ssol_i_eq_SR, p);
ssol_i_1S_SR = stabilize(ssol_i_1S_SR);

%% re-enamble figson, pinana and calcJ options

sol_eq.p.figson = 1; % re-enable figures creation for this solution
sol_eq.p.Ana = 1;
sol_eq.p.calcJ = original_p.calcJ;

ssol_eq.p.figson = 1; % re-enable figures creation for this solution
ssol_eq.p.Ana = 1;
ssol_eq.p.calcJ = original_p.calcJ;

sol_i_eq.p.figson = 1; % re-enable figures creation for this solution
sol_i_eq.p.Ana = 1;
sol_i_eq.p.calcJ = original_p.calcJ;

sol_i_eq_SR.p.figson = 1; % re-enable figures creation for this solution
sol_i_eq_SR.p.Ana = 1;
sol_i_eq_SR.p.calcJ = original_p.calcJ;

ssol_i_eq.p.figson = 1; % re-enable figures creation for this solution
ssol_i_eq.p.Ana = 1;
ssol_i_eq.p.calcJ = original_p.calcJ;

sol_eq_SR.p.figson = 1; % re-enable figures creation for this solution
sol_eq_SR.p.Ana = 1;
sol_eq_SR.p.calcJ = original_p.calcJ;

ssol_eq_SR.p.figson = 1; % re-enable figures creation for this solution
ssol_eq_SR.p.Ana = 1;
ssol_eq_SR.p.calcJ = original_p.calcJ;

ssol_i_eq_SR.p.figson = 1; % re-enable figures creation for this solution
ssol_i_eq_SR.p.Ana = 1;
ssol_i_eq_SR.p.calcJ = original_p.calcJ;

sol_1S.p.figson = 1; % re-enable figures creation for this solution
sol_1S.p.Ana = 1;
sol_1S.p.calcJ = original_p.calcJ;

sol_1S_SR.p.figson = 1; % re-enable figures creation for this solution
sol_1S_SR.p.Ana = 1;
sol_1S_SR.p.calcJ = original_p.calcJ;

sol_i_1S.p.figson = 1; % re-enable figures creation for this solution
sol_i_1S.p.Ana = 1;
sol_i_1S.p.calcJ = original_p.calcJ;

sol_i_1S_SR.p.figson = 1; % re-enable figures creation for this solution
sol_i_1S_SR.p.Ana = 1;
sol_i_1S_SR.p.calcJ = original_p.calcJ;

ssol_1S.p.figson = 1; % re-enable figures creation for this solution
ssol_1S.p.Ana = 1;
ssol_1S.p.calcJ = original_p.calcJ;

ssol_1S_SR.p.figson = 1; % re-enable figures creation for this solution
ssol_1S_SR.p.Ana = 1;
ssol_1S_SR.p.calcJ = original_p.calcJ;

ssol_i_1S.p.figson = 1; % re-enable figures creation for this solution
ssol_i_1S.p.Ana = 1;
ssol_i_1S.p.calcJ = original_p.calcJ;

ssol_i_1S_SR.p.figson = 1; % re-enable figures creation for this solution
ssol_i_1S_SR.p.Ana = 1;
ssol_i_1S_SR.p.calcJ = original_p.calcJ;

pinana(ssol_i_1S_SR); % for the last simulation draw graphics
disp('EQUILIBRATION COMPLETE')

toc

end
