function [t] = meshgen_t(p)

% v2struct(params);

% change for plotting time mesh
meshfigon = false;

% define solution mesh either logarithmically or linearly spaced points
switch p.tmesh_type 
    case 1
        t = linspace(0,p.tmax,p.tpoints);
        %xspace = linspace(0,xmax,pp+pii+pn);      % Array of point values for optical interp1

    case 2
        t = logspace(log10(p.t0),log10(p.tmax),p.tpoints) - p.t0;
        % For use with TPV only! Start with 100 points up to the decay then
        % switches to log
        
    case 3
        t = [linspace(0, p.pulsestart + p.pulselen, 0.2 * p.tpoints), p.pulsestart + p.pulselen + logspace(log10(p.pulselen/(0.8*p.tpoints*1e3)), log10(p.tmax - p.pulsestart- p.pulselen), 0.8*p.tpoints)];
        
    case 4
        t = [linspace(0, p.pulsestart, 0.1*p.tpoints), p.pulsestart+ logspace(log10((p.pulselen/(0.7*p.tpoints*1e3))), log10(p.pulselen), 0.2*p.tpoints), p.pulsestart + p.pulselen + logspace(log10(p.pulselen/(0.7*p.tpoints*1e3)), log10(p.tmax- p.pulsestart- p.pulselen), 0.8*p.tpoints)];
    otherwise % will never happen
        t = 0;
        error('tmesh_type not recognized')
end

if meshfigon

    ptmir = 1:1:length(t);
    
    figure(200);
    plot(t, ptmir, '.');
    xlabel('Time');
    ylabel('Point');

end
