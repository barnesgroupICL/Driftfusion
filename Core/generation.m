function gx = generation(par, source_type, laserlambda)
% This function calls the appropriate function (based on the chosen optical model) to calculate
% generation profiles as a function of position
% SOURCE_TYPE = either 'AM15' or 'laser'
% LASERLAMBDA = Laser wavelength - ignored if SOURCE_TYPE = AM15
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
xsolver = par.x_sub;
switch par.OM
    case 0
        % This currently results in the generation profile being stored twice and could be optimised
        gx = build_property(par.g0, xsolver, par, 'zeroed', 0);    
    case 1
        % beerlambert(par, x, source_type, laserlambda, figson)
        gx = beerlambert(par, par.xx, source_type, laserlambda, 0);
        % interpolate for i+0.5 mesh
        gx = interp1(par.xx, gx, xsolver);
end

end