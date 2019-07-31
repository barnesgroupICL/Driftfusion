function Framefile = moviemake2(solstruct, @plotfun, xrange)

% Makes a frame file F from a solution
% Currently configured to output the band diagram
% @PLOTFUN is the name of the plotting function that you wish to use
% XRANGE is a two element array with XL and XR
ionfigon = 0;
capfigon = 0;

pp1 = find(Vapp_arr < 1);
pp1 = pp1(end);

for i=1:pp1%length(Vapp_arr)
    
    fig1 = figure(30)
    clf


    F(i) = getframe(fig1);
    
    %
end


end
