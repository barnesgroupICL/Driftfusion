function varout = forward_difference(var, x)
varout = zeros(1, length(var)-1);

for i = 1:length(var)-1
    varout(1, i) = (var(i+1)-var(i))/(x(i+1)-x(i));
end

end