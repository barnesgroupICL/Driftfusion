function plot_CVs(sol_CVs, light_intensities, xlimits, ylimits)

% initialise legend
leg = [];
for i = 1:length(light_intensities)
    dfplot.JtotVapp(sol_CVs.sol(i),0)
    leg = [leg; string([num2str(light_intensities(i)), ' sun'])];
    hold on
end
legend(leg)

if xlimits ~= [0, 0]
    xlim(xlimits)
end

if ylimits ~= [0, 0] 
    ylim(ylimits)
end

hold off

end