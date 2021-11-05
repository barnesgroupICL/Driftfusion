function plot_CVs(sol_CV, light_intensity, xlimits, ylimits)

% initialise legend
leg = [];
%for i = 1:length(light_intensities)
    dfplot.JtotVapp(sol_CV, 0)
    leg = [string([num2str(light_intensity), ' sun'])];
    %hold on
%end
legend(leg)

if xlimits ~= [0, 0]
    xlim(xlimits)
end

if ylimits ~= [0, 0] 
    ylim(ylimits)
end

hold off

end