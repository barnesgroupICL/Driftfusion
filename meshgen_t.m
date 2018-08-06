function [t] = meshgen_t(p)

% define solution mesh either logarithmically or linearly spaced points
if p.tmesh_type == 1

    t = linspace(0,p.tmax,p.tpoints);
    %xspace = linspace(0,xmax,pp+pii+pn);      % Array of point values for optical interp1

elseif p.tmesh_type == 2
   
    t = logspace(log10(p.t0),log10(p.tmax),p.tpoints) - p.t0;

    % For use with TPV only! Start with 100 points up to the decay then
    % switches to log
elseif p.tmesh_type == 3
    
    t = [linspace(0, p.pulsestart+p.pulselen, 0.2*p.tpoints), p.pulsestart + p.pulselen + p.logspace(log10(p.deltat), log10(p.tmax- p.pulsestart- p.pulselen), 0.8*p.tpoints)];
    
end


if p.mesht_figon == 1

    ptmir = 1:1:length(t);
    
    figure(200);
    plot(t, ptmir, '.');
    xlabel('Position');
    ylabel('Time');

end

end