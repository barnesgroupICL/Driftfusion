function plot_sigma_light(sol_CV, light_intensities)

% initialise legend
leg = [];
for i = 1:size(sol_CVs.sigma, 1)
    Vapp = dfana.calcVapp(sol_CVs.sol(i));
    
    figure(102) 
    plot(Vapp, sol_CVs.sigma(i,:))
    hold on
    leg = [leg; string([num2str(light_intensities(i)), ' sun'])];
end
xlabel('Applied voltage [V]')
ylabel('Inferred Conductivity [S cm-1]')
legend(leg);
hold off

%% Plot sigma(Vapp=0) vs light intensity
figure(101)
plot(light_intensities, sol_CVs.sigma(:,1), '-o')
xlabel('Light intensity [suns]')
ylabel('Conductivity at Vapp = 0 V [S cm-1]')

end