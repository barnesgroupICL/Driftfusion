function gxt = gfun(gx, int, g_fun_type, par, argsin)
% Function generator for generation profile
% GX = The input generation profile as a function of position
% INT = Intensity
% G_FUN_TYPE = Type of function
%% 'constant'
% VARARGIN = INT
%% 'sweep'
% VARARGIN = {Initial_intensity, Final_intensity}
% Here INT is ignored - this was the easiest way to maintain backwards
% compatibility with different protocols and may be updated in future
% releases
%% 'square'
% VARAGIN = {low_int, high_int, time_period, duty_cycle}
% PAR = the parameters object

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
        tarr = sfg('du',t,1/argsin(3),argsin(4),0,0);
        tarr = tarr';
        gxt = tarr*gx;
    case 'sin'
        
end
