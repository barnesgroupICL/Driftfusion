function [t] = meshgen_t(params)

v2struct(params);

% define solution mesh either logarithmically or linearly spaced points
switch p.tmesh_type
    case 1
        t = linspace(0,tmax,tpoints);

    case 2
        t = logspace(log10(t0),log10(tmax),tpoints) - t0;

    % For use with TPV only! Start with 100 points up to the decay then
    % switches to log
    case 3
        t = [linspace(0, pulsestart+pulselen, 0.2*tpoints), pulsestart + pulselen + logspace(log10(pulselen/(0.8*tpoints*1e3)), log10(tmax- pulsestart- pulselen), 0.8*tpoints)];

    case 4
        t = [linspace(0, pulsestart, 0.1*tpoints), pulsestart+ logspace(log10((pulselen/(0.7*tpoints*1e3))), log10(pulselen), 0.2*tpoints), pulsestart + pulselen + logspace(log10(pulselen/(0.7*tpoints*1e3)), log10(tmax- pulsestart- pulselen), 0.8*tpoints)];
    otherwise
        error('DrIFtFUSION:tmesh_type', [mfilename ' - tmesh_type not recognized'])
end

end
