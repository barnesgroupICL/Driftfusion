function anasdp(sdpsol, Jtr_time)
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
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