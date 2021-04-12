function EA_script_ana(type, filename, varargin)
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
            case '1h'
                plot(varargin{isol}.Freq, varargin{isol}.AC_ExDC_E_amp, options{:})
            case '2h'
                plot(varargin{isol}.Freq, varargin{isol}.AC_Efield2_amp, options{:})
            case 'phase'
                phase_n_deg = rad2deg(wrapTo2Pi(varargin{isol}.AC_ExDC_E_phase));
                plot(varargin{isol}.Freq, phase_n_deg, options{:})
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
    
    ax = gca;
    ax.XScale = 'log'; % for putting the scale in log
    ax.YScale = 'log'; % for putting the scale in log
    xlabel('Frequency [Hz]');
    ylabel('Abs(E_{AC}^2) [V^2/cm^2]');

legend(legendarr, 'Location', 'southeast')
legend boxoff
if 7~=exist(filename,'dir')
    mkdir(filename)
end
saveas(fig, [char(filename) filesep char(filename) char(['-EA_' type '.fig'])])
saveas(fig, [char(filename) filesep char(filename) char(['-EA_' type '.png'])])

end