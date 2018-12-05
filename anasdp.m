function anasdp(sdpsol, Jtr_time)

% sdpsol.t_Jtr
% sdpsol.Jtr
% sdpsol.tdwell_arr
p1 = find(sdpsol.t_Jtr >= Jtr_time);
p1 = p1(1);

sdpsol.Jtr_time_arr = sdpsol.Jtr(p1,:);

figure(1)
semilogx(sdpsol.tdwell_arr, sdpsol.Jtr_time_arr)
xlabel('t_{dwell} [s]')
ylabel('Jtr [Acm-2]')

figure(2)

for i = 1:length(sdpsol.tdwell_arr)
    
    plot(sdpsol.t_Jtr, sdpsol.Jtr(:,i))
    xlabel('time [s]')
    ylabel('J [Acm-2]')
    hold on    
end
hold off

end