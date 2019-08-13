function gxt = fun_gen(gx, A, fun_type, argsin, t, tmax)
% Function generator for generation profile
% GX = The input profile as a function of position for 2 dimensional
% functions
% A = Amplitude
% fun_type = Type of function
% PAR = the parameters object
%% 'constant'
% ARGSIN is ignored and A is used
%% 'sweep'
% ARGSIN = {A_initial, A}
% Here A is ignored - this was the easiest way to maintain backwards
% compatibility with different protocols and may be updated in future
% releases
%% 'square'
% ARGSIN = {A_low, A_high, time_period, duty_cycle}
%% 'sin'
% ARGSIN = {DC_int, Delta_int, frequency}

switch fun_type
    case 'constant'
        gxt = A*repmat(gx, length(t), 1);
    case 'sweep'
        range = argsin;
        for i=1:length(t)
            gxt(i,:) = gx*(range(1) + ((range(2)-range(1))*t(i)*(1/tmax)));
        end
    case 'square'
        % Generate intensity array
        A_arr = argsin(1) + (argsin(2)-argsin(1))*sfg('du',t,1/argsin(3),argsin(4),0,0);
        A_arr = A_arr';
        gxt = A_arr*gx;
    case 'sin'
        A_arr = argsin(1) + argsin(2)*(2*sin(2*pi*argsin(3)*t));
        A_arr = A_arr';
        gxt = A_arr*gx;
end
