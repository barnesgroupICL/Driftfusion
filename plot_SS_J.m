dfplot.Jtott(sol_relax.intrinsic_el,0)
hold on
dfplot.Jtott(sol_relax.intrinsic_ion,0)
hold on
dfplot.Jtott(sol_relax.hom_el,0)
hold on
dfplot.Jtott(sol_relax.hom_ion,0)
hold off
ylim([0.0216, 0.022])
legend('Intrinsic', 'Intrinsic with Shottky Defects', 'Homojunction', 'Homojunction with Shottky Defects')