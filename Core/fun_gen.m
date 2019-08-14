function fun = fun_gen2(g_fun_type)
% Function generator for generation profile
% GX = The input profile as a function of position for 2 dimensional
% functions
% A = Amplitude
% G_FUN_TYPE = Type of function
% PAR = the parameters object
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
% COEFF = [DC_int, Delta_A, frequency, phase]

switch g_fun_type
    case 'constant'
        fun = @(coeff, t) coeff(1);
    case 'sweep'
        fun = @(coeff, t) coeff(1) + (coeff(2)-coeff(1))*t/coeff(3);
    case 'square'
        % Generate intensity array
        fun = @(coeff, t) coeff(1) + (coeff(2)-coeff(1))*lt(mod(t,coeff(3))*1/coeff(3),coeff(4)/100);
    case 'sin'
        fun = @(coeff, t) coeff(1) + coeff(2)*(2*sin(2*pi*coeff(3)*t + coeff(4)));
end
