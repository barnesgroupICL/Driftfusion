function varihalf = getvarihalf(var)

for i = 1:length(var)-1
    varihalf(i) = var(i)+((var(i+1)-var(i))/2);
end
end