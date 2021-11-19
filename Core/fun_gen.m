function fun = fun_gen(fun_type)
% Function generator for generation profile
% Some of the functions in this code are based on SFG.M by Hiroyuki Kato Copyright (c) 2013
% The license details are contained with sfg_license
% G_FUN_TYPE = Type of function
%% 'constant'
% COEFF = [Amplitude]
%% 'sweep'
% COEFF = [Amplitude_initial, Amplitude_final, tmax]
% Here A is ignored - this was the easiest way to maintain backwards
% compatibility with different protocols and may be updated in future
% releases
%% 'square'
% COEFF = [A_low, A_high, time_period, duty_cycle]
%% 'sin'
% COEFF = [DC_offset, Delta_AC, frequency, phase]
%% 'sweepAndStill'
% COEFF = [A_start, A_end, sweep_duration]
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Start code
switch fun_type
    case 'constant'
        fun = @(coeff, t) coeff(1);
    case 'sweep'
        fun = @(coeff, t) coeff(1) + (coeff(2)-coeff(1))*t/coeff(3);
    case 'square'
        % Generate intensity array
        fun = @(coeff, t) coeff(1) + (coeff(2)-coeff(1))*lt(mod(t,coeff(3))*1/coeff(3),coeff(4)/100);
    case 'sin'
        fun = @(coeff, t) coeff(1) + coeff(2)*(sin(2*pi*coeff(3)*t + coeff(4)));
    case 'sweep2'
        fun = @(coeff, t) coeff(1) + [(coeff(2)-coeff(1))*t(t<=coeff(3)/2)/coeff(3),...
                                   -(coeff(2)-coeff(1))*(t(t>coeff(3)/2)/coeff(3) -  2*t(t==coeff(3)/2)/coeff(3))];
    case 'tri'
        % COEFF = [OFFSET, V1, V2, periods, tperiod]  tmax is defined by the input
        fun = @(coeff, t) triangle_fun(coeff, t);
    case 'sweepAndStill'
        % COEFF = [A_start, A_end, sweep_duration]
        % do a sweep from coeff(1) to coeff(2) with the duration of
        % coeff(3), then stay stable at coeff(2)
        fun = @(coeff, t) coeff(2) + (coeff(1)-coeff(2))*max(0,1-t/coeff(3));
end
