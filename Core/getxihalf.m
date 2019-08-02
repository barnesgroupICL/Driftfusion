function xsolver = getxihalf(sol)
x = sol.par.xx;
for i = 1:length(x)-1
    xsolver(i) = x(i)+((x(i+1)-x(i))/2);
end
end