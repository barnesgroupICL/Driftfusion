function Vout = triangle_fun_singlecycle(coeff, t)
% Generates a cyclic voltage function

V0 = coeff(1);
V1 = coeff(2);
V2 = coeff(3);
tmax = coeff(4);

Vout = zeros(1, length(t));
% scan rate
deltaV = abs(V1-V0)+abs(V2-V1)+abs(V0-V2);
k = deltaV/tmax;

t1 = abs(V1-V0)/k;
t2 = t1+abs(V2-V1)/k;
t3 = t2+abs(V0-V2)/k;

Vout(t<t1) = V0+(V1-V0).*t(t<t1)./t1;
Vout(t>=t1 & t<t2) = V1+((V2-V1).*(t(t>=t1 & t<t2)-t1))./(t2-t1);
Vout(t>=t2) = V2+(V0-V2).*(t(t>=t2)-t2)/(tmax-t2);

end