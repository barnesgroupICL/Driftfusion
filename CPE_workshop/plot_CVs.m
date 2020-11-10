function plot_CVs(sol_CVs, light_intensities, xlimits, ylimits)

% initialise legend
leg = [];
for i=1:length(light_intensities)
    dfplot.JtotVapp(sol_CVs(i),0)
    leg = [leg; string([num2str(light_intensities(i)), ' sun'])];
    hold on
end
legend(leg)
xlim(xlimits)
ylim(ylimits)
hold off

end