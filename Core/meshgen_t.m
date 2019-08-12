function [t] = meshgen_t(par)
%% Creates the output time mesh
% Note- ODE15S uses an adaptive time step and then interpolates to the
% requested times

switch par.tmesh_type
    case 1
        t = linspace(0,par.tmax,par.tpoints);     
    case 2
        t = logspace(log10(par.t0),log10(par.tmax),par.tpoints);
        %To avoid rounding errors
        t(1) = 0;
        t(end) = par.tmax;
    case 3
        % For use with TPV only! Start with 100 points up to the decay then
        % switches to log
        t = [linspace(0, par.pulsestart+par.pulselen, 0.2*par.tpoints), par.pulsestart + par.pulselen + par.logspace(log10(par.deltat), log10(par.tmax- par.pulsestart- par.pulselen), 0.8*par.tpoints)];
    case 4
        % 2 log meshes consecutively
        t1 = logspace(log10(par.t0),log10(par.tmax/2),round(par.tpoints/2));
        t1(1) = 0;
        t1(end) = par.tmax/2;
        t2 = t1(end) + logspace(log10(par.t0),log10(par.tmax/2),round(par.tpoints/2));
        t2(end) = par.tmax;
        t = [t1, t2];
end

if par.mesht_figon == 1
    
    ptmir = 1:1:length(t);
    
    figure(200);
    plot(t, ptmir, '.');
    xlabel('Position');
    ylabel('Time'); 
end

end