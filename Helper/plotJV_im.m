function plotJV_im(sol)

plot(sol.V, sol.J*1e-3)
xlabel('Applied Voltage [V]')
ylabel('Current Density [mA cm-2]')

end