function IS_script_ana(type, filename, varargin)
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% sdpsol.t_Jtr
% sdpsol.Jtr
% sdpsol.tdwell_arr

fig = figure('Name', 'Amplitude of EA second harmonic E_{AC}^2', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
    hold off
    %ymin = Inf;
    %ymax = -Inf;
    for i = 1:(length(varargin)/3)
        isol = i*3-2;
        legendarr{i} = varargin{isol+1};
        options = varargin{isol+2};
        switch type
            case 'impedance_re'
                h(i) = plot(varargin{isol}.Freq, varargin{isol}.impedance_re, options{:});
                xlabel('Frequency [Hz]');
                ylabel('Re(Z) [\Omega cm^2]');
                ax = gca;
                ax.XScale = 'log'; % for putting the scale in log
                ax.YScale = 'log'; % for putting the scale in log
            case 'impedance_im'
                h(i) = plot(varargin{isol}.Freq, varargin{isol}.impedance_im, options{:});
                xlabel('Frequency [Hz]');
                ylabel('-Im(Z) [\Omega cm^2]');
                ax = gca;
                ax.XScale = 'log'; % for putting the scale in log
                ax.YScale = 'log'; % for putting the scale in log
            case 'impedance_abs'
                h(i) = plot(varargin{isol}.Freq, varargin{isol}.impedance_abs, options{:});
                xlabel('Frequency [Hz]');
                ylabel('Abs(Z) [\Omega cm^2]');
                ax = gca;
                ax.XScale = 'log'; % for putting the scale in log
                ax.YScale = 'log'; % for putting the scale in log
            case 'capacitance'
                h(i) = plot(varargin{isol}.Freq, varargin{isol}.cap, options{:});
                hold on
                plot(varargin{isol}.Freq, -varargin{isol}.cap, options{:}, 'Marker', 'o', 'MarkerSize', 7)
                xlabel('Frequency [Hz]');
                ylabel('\omega^{-1} x Im(Z^{-1}) [F/cm^2]');
                ax = gca;
                ax.XScale = 'log'; % for putting the scale in log
                ax.YScale = 'log'; % for putting the scale in log
            case 'phase'
                phase_n_deg = rad2deg(varargin{isol}.Jtot_phase);
                h(i) = plot(varargin{isol}.Freq, -phase_n_deg, options{:});
                xlabel('Frequency [Hz]');
                ylabel('Z Phase [deg]');
                ax = gca;
                ax.XScale = 'log'; % for putting the scale in log
            case 'nyquist'
                h(i) = plot(varargin{isol}.impedance_re, -varargin{isol}.impedance_im, options{:});
                xlabel('Re(Z) [\Omega cm^2]');
                ylabel('-Im(Z) [\Omega cm^2]');
            otherwise
                error([mfilename(1) ' - *type* option needed'])
        end
        %ymin = min(ymin, min());
        %ymax = max(ymax, max());
        %range = ymax-ymin;
        %ylim([ymin-0.03*range, ymax+0.03*range])
        %xlim([min(min(EA_results.Freq)), max(max(EA_results.Freq))])

        hold on
    end
    

legend(h, legendarr, 'Location', 'northeast')
legend boxoff
if 7~=exist(filename,'dir')
    mkdir(filename)
end
saveas(fig, [char(filename) filesep char(filename) char(['-IS_' type '.fig'])])
saveas(fig, [char(filename) filesep char(filename) char(['-IS_' type '.png'])])

end