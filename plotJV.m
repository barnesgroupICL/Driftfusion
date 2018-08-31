function plotJV(JV, option)

% option        1 = dark only, 2 = light only, 3 = dark & light
% JV is a structure containing dark and illuminated JVs

figure(11)
ylim([-30, 30]);
hold on

if option == any([1, 3])
    pinana(JV.dk.f);
    pinana(JV.dk.r);
    
end

if option == any([2, 3])
    pinana(JV.ill.f);
    pinana(JV.ill.r);
end

hold off

end