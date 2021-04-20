function SDP_script_ana(Jtr_time, filename, varargin)
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

% 
% SDP_script_ana(1e-7, '20210302_spiro_CBTiO2',...
%     sdpsol_spiro_lowerCBTiO2_alt_dark, 'lower CBTiO2 dark',{':k','LineWidth',3},...
%     sdpsol_spiro_lowerCBTiO2_alt_01sun, 'lower CBTiO2 0.1 sun',{':b','LineWidth',3},...
%     sdpsol_spiro_lowerCBTiO2_alt_1sun, 'lower CBTiO2 1 sun',{':r','LineWidth',3},...
%     sdpsol_spiro_alt_dark, ' dark',{'--k'},...
%     sdpsol_spiro_alt_01sun, ' 0.1 sun',{'--b'},...
%     sdpsol_spiro_alt_1sun, ' 1 sun',{'--r'},...
%     sdpsol_spiro_higherCBTiO2_alt_dark, 'higher CBTiO2 dark',{'-k'},...
%     sdpsol_spiro_higherCBTiO2_alt_01sun, 'higher CBTiO2 0.1 sun',{'-b'},...
%     sdpsol_spiro_higherCBTiO2_alt_1sun, 'higher CBTiO2 1 sun',{'-r'})



for i = 1:(length(varargin)/3)
    isol = round((i*3)-2);
    p1_list = find(varargin{isol}.t_Jtr >= Jtr_time);
    p1(i) = p1_list(1);
    Jtr_time_arr(i,:) = varargin{isol}.Jtr(p1(i),:);%+varargin{isol}.Jdk(p1(i));
    if ~ismembertol(i, 1)
        if ~ismembertol(varargin{isol-3}.t_Jtr(p1(i-1)), varargin{isol}.t_Jtr(p1(i)))
            warning('Solutions are plotted at different times!!!!')
            disp(varargin{isol-3}.t_Jtr(p1(i-1)))
            disp(varargin{isol}.t_Jtr(p1(i)))
        end
    end
end

fig = figure('Name', 'SDP vs dwell time', 'NumberTitle', 'off', 'units','normalized', 'outerposition',[0 0 1 1]);
hold off
ymin = Inf;
ymax = -Inf;
%legendarr = zeros(length(varargin)/3, 1);
for i = 1:(length(varargin)/3)
    isol = i*3-2;
    legendarr{i} = varargin{isol+1};
    options = varargin{isol+2};
    plot(varargin{isol}.tdwell_arr, Jtr_time_arr(i,:), options{:})
    ymin = min(ymin, min(Jtr_time_arr(i,(end/2):end)));
    ymax = max(ymax, max(Jtr_time_arr(i,1:(end/2))));
    hold on
end
xlabel('t_{dwell} [s]')
ylabel('Jtr [Acm-2]')
ax = gca;
ax.XScale = 'log'; % for putting the scale in log
%xlim([2e-8, 1e3])
range = ymax-ymin;
ylim([ymin-0.03*range, ymax+0.03*range])
legend(legendarr)
legend boxoff
if 7~=exist(filename,'dir')
    mkdir(filename)
end
saveas(fig, [char(filename) filesep char(filename) char('-SDP.fig')])
saveas(fig, [char(filename) filesep char(filename) char('-SDP.png')])

end
