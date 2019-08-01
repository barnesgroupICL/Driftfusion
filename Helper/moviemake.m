function Framefile = moviemake2(sol, plotfun, xrange, yrange, movie_name, Vcounter, tcounter)
% Makes a frame file F from a solution
% Currently configured to output the band diagram
% @PLOTFUN is the name of the plotting function that you wish to use
% XRANGE is a two element array with XL and XR in cm
% yRANGE is a two element array with YLOW and YHIGH
% MOVIE_NAME is the desired out name


ionfigon = 0;
capfigon = 0;

Vapp = dfana.calcVapp(sol);

for i = 1:length(sol.t)
    %figure(600)
    clf
    
    plotfun(sol, sol.t(i))
    fig1 = gca;
    subplot(2,1,1);
    ylim([[-0.7, 0.1]])
    
    subplot(2,1,2);
    ylim([0e18, 2e18])

%     subplot(3,1,3);
%     ylim([0,2e18])
    
    if xrange ~= 0
        xlim([xrange(1)*1e7, xrange(2)*1e7])
    end
    
    if yrange ~= 0
        ylim([yrange(1), yrange(2)])
    end
    
    Framefile(i) = getframe(gcf);
    
    % Voltage counter
    if Vcounter
        dim = [.2 0 .3 .3];
        Vnow = round(Vapp(i), 2, 'decimal');
        anno = ['V = ', num2str(Vnow), ' V'];
        T = annotation('textbox', dim, 'String', anno, 'FitBoxToText','on');
        T.FontSize = 16;
        drawnow
    end
    
    % Time counter
    if tcounter
        dim = [.2 0 .3 .3];
        tnow = round(sol.t(i), 2, 'decimal');
        anno = ['t = ', num2str(tnow), ' s'];
        T = annotation('textbox', dim, 'String', anno, 'FitBoxToText','on');
        T.FontSize = 16;
        drawnow
    end
    
end

moviewrite(Framefile, movie_name);

function moviewrite(Framefile, movie_name)
% Write a frame file to an avi movie
% name is a string with the desired filename- do NOT include .avi extension

% name is a string for the final video
movie_name = [movie_name, '.avi'];

% Write to file
myVideo = VideoWriter(movie_name);
myVideo.FrameRate = 6;
open(myVideo);
writeVideo(myVideo, Framefile);
close(myVideo);

end

end
