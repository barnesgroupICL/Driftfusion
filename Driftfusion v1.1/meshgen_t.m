function [t] = meshgen_t(p)

% define solution mesh either logarithmically or linearly spaced points
switch p.tmesh_type
    case 1
        t = linspace(0, p.tmax, p.tpoints);

    case 2
        t = logspace(log10(p.t0), log10(p.tmax), p.tpoints) - p.t0;

    % For use with TPV only! Start with 100 points up to the decay then
    % switches to log
    case 3
        t = [linspace(0, p.pulsestart+p.pulselen, 0.2*p.tpoints),...
            p.pulsestart + p.pulselen + logspace(log10(p.pulselen/(0.8*p.tpoints*1e3)), log10(p.tmax-p.pulsestart-p.pulselen), 0.8*p.tpoints)];

    case 4
        t = [linspace(0, p.pulsestart, 0.1*p.tpoints),...
            p.pulsestart + logspace(log10(p.pulselen/(0.7*p.tpoints*1e3)), log10(p.pulselen), 0.2*p.tpoints),...
            p.pulsestart + p.pulselen + logspace(log10(p.pulselen/(0.7*p.tpoints*1e3)), log10(p.tmax-p.pulsestart-p.pulselen), 0.8*p.tpoints)];
    otherwise
        error('DrIFtFUSION:tmesh_type', [mfilename ' - tmesh_type not recognized'])
end

end
