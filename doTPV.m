function [ssol_TPV] = doTPV(ssol_eq, mui, eqtime, biasint, pulseint, pulselen, decaytime)

% This function uses the open circuit equilibrium solution to perform a
% transient photovoltage measure
%% Input arguments
% ssol_eq   = OC equilibrium solution
% mui       = ion mobility [cm2V-1s-1]
% eqtime    = equilibration time [s]
% biasint   = bias light intensity [Suns eq.]
% pulseint  = pulse light intensity [Sun eq.]
% pulselen  = Pulse length [s]
% decaytime = Length of recorded decay [s]

% Read-in parameters
p = ssol_eq.params;

%% Initial parameters

p.OC = 1;                 % Closed circuit = 0, Open Circuit = 1 
p.Int = biasint;          % Bias Light intensity (Suns Eq.)
p.tmax = 1e-4;            % Time
p.t0 = p.tmax/1e8;
p.tmesh_type = 2;       
p.pulseon = 0;            % Switch pulse on TPC or TPV
p.Vapp = 0;               % Applied bias
p.figson = 1;             % Toggle figures on/off
p.meshx_figon = 0;        % Toggles x-mesh figures on/off
p.mesht_figon = 0;        % Toggles t-mesh figures on/off
p.calcJ = 0;              % Calculates Currents- slows down solving calcJ = 1, calculates DD currents at every position, calcJ = 2, calculates DD at boundary.
p.JV = 0;                 % Toggle run JV scan on/off
p.Ana = 1;                % Toggle on/off analysis

%% Initial solution under illumination
disp('Solving for device under illumination..')
ssol_Il = pindrift(ssol_eq, p);

disp('Complete.')

disp('Equilibration...')
%% Switch on ion mobility
p.mui = mui;
p.tmax = eqtime;
p.t0 = eqtime/1e6;

ssol_i_Il = pindrift(ssol_Il, p);
disp('Complete.')

%% Do the TPV
disp('Doing the TPV')
p.pulseon = 1;
p.pulselen = pulselen;       % Transient pulse length
p.pulseint = pulseint;          % Transient pulse intensity- for BL and TM models, 100 mW/cm2 assumed
p.pulsestart = 1e-7;     % Time before pulse- required for baseline
p.tmesh_type = 3;
p.tmax = decaytime;
p.t0 = p.tmax/1e6;
p.tpoints = 400;
p.deltat = p.tmax/(1e6*p.tpoints);

ssol_TPV = pindrift(ssol_i_Il, p);
disp('Complete.')

disp('doTPV complete.')

end
