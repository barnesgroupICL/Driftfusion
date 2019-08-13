function gxt = fun_gen(gx, A, g_fun_type, par, argsin)
% Function generator for generation profile
% GX = The input profile as a function of position for 2 dimensional
% functions
% A = Amplitude
% G_FUN_TYPE = Type of function
% PAR = the parameters object
%% 'constant'
% ARGSIN is ignored and INT is used
%% 'sweep'
% ARGSIN = {Initial_intensity, Final_intensity}
% Here INT is ignored - this was the easiest way to maintain backwards
% compatibility with different protocols and may be updated in future
% releases
%% 'square'
% ARGSIN = {A_low, A_high, time_period, duty_cycle}
%% 'sin'
% ARGSIN = {DC_int, Delta_int, frequency}

t = meshgen_t(par);
switch g_fun_type
    case 'constant'
        gxt = int*repmat(gx, length(t), 1);
    case 'sweep'
        range = argsin;
        for i=1:length(t)
            gxt(i,:) = gx*(range(1) + ((range(2)-range(1))*t(i)*(1/t(end))));
        end
    case 'square'
        % Generate intensity array
        int_arr = argsin(1) + (argsin(2)-argsin(1))*sfg('du',t,1/argsin(3),argsin(4),0,0);
        int_arr(end) = argsin(1);
        int_arr = int_arr';
        gxt = int_arr*gx;
    case 'sin'
        int_arr = argsin(1) + argsin(2)*(2*sin(2*pi*argsin(3)*t));
        int_arr(1) = argsin(1);     %Ensure starts at int_intial
        int_arr = int_arr';
        gxt = int_arr*gx;
end
