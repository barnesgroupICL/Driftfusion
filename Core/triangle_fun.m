function Vout = triangle_fun(coeff, t)
% Generates a triangle wave form

V0 = coeff(1);
V1 = coeff(2);
V2 = coeff(3);
cycles = coeff(4);
tmax = coeff(5);
if cycles == 1
    p_per_cycle = length(t)-1;
else
    p_per_cycle = floor(length(t)/cycles);   % Points per cycle
end

Vout = zeros(1, length(t));
% scan rate
deltaV = abs(V1-V0)+abs(V2-V1)+abs(V0-V2);
k = (deltaV)/(tmax/cycles);

t1 = abs(V1-V0)/k;
t2 = t1+abs(V2-V1)/k;
t3 = t2+abs(V0-V2)/k;

V_1cycle(t<t1) = V0+(V1-V0).*t(t<t1)./t1;
V_1cycle(t>=t1 & t<t2) = V1+((V2-V1).*(t(t>=t1 & t<t2)-t1))./(t2-t1);
V_1cycle(t>=t2) = V2+(V0-V2).*(t(t>=t2)-t2)/(tmax/cycles-t2);

for i=1:cycles
    p1 = 1+((i-1)*p_per_cycle);
    p2 = 1+(i*p_per_cycle);
    Vout(p1:p2) = V_1cycle(1:1+p_per_cycle);
end

end