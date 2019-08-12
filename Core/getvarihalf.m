function varihalf = getvarihalf(var)

varihalf = zeros(1, length(var)-1);
for i = 1:length(var)-1
    varihalf(1,i) = var(1,i)+((var(1,i+1)-var(1,i))/2);
end
end