function [t] = meshgen_t(params)

v2struct(params);

meshfigon = 0;

% define solution mesh either logarithmically or linearly spaced points
if tmesh_type == 1

    t = linspace(0,tmax,tpoints);
    %xspace = linspace(0,xmax,pp+pii+pn);      % Array of point values for optical interp1

elseif tmesh_type == 2
   
    t = logspace(log10(t0),log10(tmax),tpoints) - t0;

    % For use with TPV only! Start with 100 points up to the decay then
    % switches to log
elseif tmesh_type == 3
    
    t = [linspace(0, pulsestart+pulselen, 0.2*tpoints), pulsestart + pulselen + logspace(log10(pulselen/(0.8*tpoints*1e3)), log10(tmax- pulsestart- pulselen), 0.8*tpoints)];
   
elseif tmesh_type == 4
    
    t = [linspace(0, pulsestart, 0.1*tpoints), pulsestart+ logspace(log10((pulselen/(0.7*tpoints*1e3))), log10(pulselen), 0.2*tpoints), pulsestart + pulselen + logspace(log10(pulselen/(0.7*tpoints*1e3)), log10(tmax- pulsestart- pulselen), 0.8*tpoints)];
    
end

pt = length(t);

if meshfigon == 1

    tmir = t;
    ptmir = 1:1:length(t);
    
    figure(200);
    plot(tmir, ptmir, '.');
    xlabel('Time');
    ylabel('Point');

end

end