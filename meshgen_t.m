function [t] = meshgen_t(par)

% define solution mesh either logarithmically or linearly spaced points
if par.tmesh_type == 1

    t = linspace(0,par.tmax,par.tpoints);
    %xspace = linspace(0,xmax,pp+pii+pn);      % Array of point values for optical interp1

elseif par.tmesh_type == 2
   
    t = logspace(log10(par.t0),log10(par.tmax),par.tpoints) - par.t0;

    % For use with TPV only! Start with 100 points up to the decay then
    % switches to log
elseif par.tmesh_type == 3
    
    t = [linspace(0, par.pulsestart+par.pulselen, 0.2*par.tpoints), par.pulsestart + par.pulselen + par.logspace(log10(par.deltat), log10(par.tmax- par.pulsestart- par.pulselen), 0.8*par.tpoints)];
    
end


if par.mesht_figon == 1

    ptmir = 1:1:length(t);
    
    figure(200);
    plot(t, ptmir, '.');
    xlabel('Position');
    ylabel('Time');

end

end