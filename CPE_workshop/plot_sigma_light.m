function plot_sigma_light(sol_CVs, light_intensities, sigma)

% plot sigma vs Vapp

% initialise legend
leg = [];
for i = 1:size(sigma, 1)
    Vapp = dfana.calcVapp(sol_CVs(i));
    
    figure(102) 
    plot(Vapp, sigma(i,:))
    hold on
    leg = [leg; string([num2str(light_intensities(i)), ' sun'])];
end
xlabel('Applied voltage [V]')
ylabel('Inferred Conductivity [S cm-1]')
legend(leg);
hold off

%% Plot sigma(Vapp=0) vs light intensity
figure(101)
plot(light_intensities, sigma(:,1), '-o')
xlabel('Light intensity [suns]')
ylabel('Conductivity at Vapp = 0 V [S cm-1]')

end