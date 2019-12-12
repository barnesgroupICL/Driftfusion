function coeff = DA_sin(par, V_fun_arg, periods, ppperiod, tarr, figson)
% PPPEROID = Point per period

par.V_fun_type = 'sin';
par.V_fun_arg = V_fun_arg;
Vac = V_fun_arg(2);
Freq = V_fun_arg(3);

% Setup time mesh
par.tmesh_type = 1;
par.tmax = periods/Freq;
par.t0 = 0;
par.tpoints = ppperiod*periods;

[tau, Cap, rhomat, Fmat, Phimat, Q, J, tout] = depletion_approx_modelX1_fungen(par, tarr, figson);

% tfit = tout(round(par.tpoints/periods):end);
% tfit = tfit-tfit(1);
% Jfit = J(round(par.tpoints/periods):end);

coeff = ISwave_EA_single_demodulation(tout, J', par);

Zmag = Vac/coeff(2);
Cap = sin(coeff(3))/(2*pi*10*Zmag);

end