function [JV_1S_f, JV_1S_r] = tempscan

Tarr = [200, 250, 300, 350];
p = pc_varT;

for i = 1:length(Tarr)
    
    p.T = Tarr(i)
    [sol_eq, sol_i_eq, sol_i_eq_SR] = equilibrate(p);
   
    [~, ~, JV_1S_f(i), JV_1S_r(i)] = doJV(sol_i_eq_SR, 10, 100, 0, 0, 1.2, 2);

    figure(500)
    plot(JV_1S_f(i).Vapp', JV_1S_f(i).Jtotr(:, end), JV_1S_r(i).Vapp', JV_1S_r(i).Jtotr(:, end))
    xlabel('V [V]')
    ylabel('J [mAcm-2]')
    hold on
    
end

hold off

end
