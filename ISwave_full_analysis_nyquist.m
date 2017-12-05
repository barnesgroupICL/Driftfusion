function ISwave_full_analysis_nyquist(ISwave_struct)
% plot the Nyquist graph for Impedance Spectroscopy (IS) in a range of background light intensities,
% capacitance versus voltage oscillation frequency at various light bias

Int_colors = colormap(jet(length(ISwave_struct.Int)));
legend_Int = string(ISwave_struct.Int);

figure('Name', 'Nyquist plot of IS at various light intensities', 'NumberTitle', 'off')
    hold off
    for i = 1:length(ISwave_struct.Int)
        plot(ISwave_struct.impedance_re(i, :), ISwave_struct.impedance_im(i, :)', 'Color', Int_colors(i, :), 'Marker', 'd');
        hold on
    end
    xlabel('Re(Z) [\Omega cm^2]');
    ylabel('Im(Z) [\Omega cm^2]');
    legend(legend_Int)
    
