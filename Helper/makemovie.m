function Framefile = makemovie(sol, plotfun, xrange, yrange, movie_name, Vcounter, tcounter)
% Makes a frame file F from a solution
% Currently configured to output the band diagram
% @PLOTFUN is the name of the plotting function that you wish to use
% XRANGE is a two element array with XL and XR in cm - enter 0 for auto
% yRANGE is a two element array with YLOW and YHIGH - enter 0 for auto
% MOVIE_NAME is the desired out name
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Start code
Vapp = dfana.calcVapp(sol);

for i = 1:length(sol.t)
    clf
    plotfun(sol, sol.t(i))
    fig1 = gca;

%% You can include limits for subplots here
%     subplot(2,1,1);
%     ylim([-2, 0.2])
%
%     subplot(2,1,2);
%     ylim([0, 3.5e16])

    if xrange ~= 0
        xlim([xrange(1)*1e7, xrange(2)*1e7])
    end

    if yrange ~= 0
        ylim([yrange(1), yrange(2)])
    end
    set(gcf,'color','w');

    % Voltage counter
    if Vcounter
        dim = [.2 0 .3 .2];
        Vnow = round(Vapp(i), 2, 'decimal');
        anno = ['V = ', num2str(Vnow), ' V'];
        T = annotation('textbox', dim, 'String', anno, 'FitBoxToText','on');
        T.EdgeColor = 'none';
        T.FontSize = 16;
        drawnow
    end

    % Time counter
    if tcounter
        dim = [.2 0 .3 .3];
        tnow = round(sol.t(i), 2, 'decimal');
        anno = ['t = ', num2str(tnow), ' s'];
        T = annotation('textbox', dim, 'String', anno, 'FitBoxToText','on');
        T.EdgeColor = 'none';
        T.FontSize = 16;
        drawnow
    end

    Framefile(i) = getframe(gcf);
end

moviewrite(Framefile, movie_name);

function moviewrite(Framefile, movie_name)
% Write a frame file to an avi movie
% name is a string with the desired filename- do NOT include .avi extension

% name is a string for the final video
movie_name = [movie_name, '.avi'];

% Write to file
myVideo = VideoWriter(movie_name);
myVideo.FrameRate = 20;
open(myVideo);
writeVideo(myVideo, Framefile);
close(myVideo);

end

end
