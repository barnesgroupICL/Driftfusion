function plotJV(JV, option)
% JV - a soultion from doJV
% OPTION - 1 = dark only, 2 = light only, 3 = dark & light
% JV is a structure containing dark and illuminated JVs

if option == 1 || option == 3
    JV.dk.f.par.JV = 1;
    JV.dk.r.par.JV = 1;
    
    dfana(JV.dk.f);
    figure(11)
    %xlim([-0.2, 1.15])
    ylim([-25e-3, 5e-3]);
    hold on
    dfana(JV.dk.r);

    JV.dk.f.par.JV = 0;
    JV.dk.r.par.JV = 0;
end

if option == 2 || option == 3
    JV.ill.f.par.JV = 1;
    JV.ill.r.par.JV = 1;
    
    dfana(JV.ill.f);
    
    figure(11)
    %xlim([-0.2, 1.15])
    ylim([-25e-3, 5e-3]);
    hold on
    dfana(JV.ill.r);
    
    JV.ill.f.par.JV = 0;
    JV.ill.r.par.JV = 0;
end

hold off

end