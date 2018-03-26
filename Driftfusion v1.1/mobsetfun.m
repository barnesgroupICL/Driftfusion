function [p] = mobsetfun(mob, mobi, p)

% function to switch on and off mobilities. Execute like this:
% p = mobsetfun(1, params)

% params = pinParams

    p.mue_p = mob;
    p.muh_p = mob;
    p.mue_i = mob;
    p.muh_i = mob;
    p.mue_n = mob;
    p.muh_n = mob;
    
    p.mui = mobi;

end