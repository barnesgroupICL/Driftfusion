function plotJV(JV, option)

% option        1 = dark only, 2 = light only, 3 = dark & light
% JV is a structure containing dark and illuminated JVs

figure(11)
ylim([-25, 5]);
hold on

if option == 1 || option == 3
    pinana(JV.dk.f);
    pinana(JV.dk.r);
    
end

if option == 2 || option == 3
    pinana(JV.ill.f);   
    pinana(JV.ill.r);
end

hold off

end