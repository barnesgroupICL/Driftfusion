function y = triangle_fun(coeff, t)
% Generates a triangle wave form
% Readout coefficients:
y0 = coeff(1);
y1 = coeff(2);
y2 = coeff(3);
periods = coeff(4);
tperiod = coeff(5);
tmax = tperiod*periods;

deltay = abs(y1-y0)+abs(y2-y1)+abs(y0-y2);
k = deltay/tperiod;

t0 = 0;
t1 = abs(y1-y0)/k;
t2 = t1+abs(y2-y1)/k;

y = lt(mod(t,tperiod), t1).*(y0+(y1-y0).*mod(t,tperiod)./t1) +...
    ge(mod(t,tperiod), t1).*lt(mod(t,tperiod), t2).*(y1+((y2-y1).*mod((t-t1),tperiod)./(t2-t1))) +...
    ge(mod(t,tperiod), t2).*lt(mod(t,tperiod), tmax).*(y2+((y0-y2).*mod((t-t2),tperiod)./(tperiod-t2)));

end