function gx = generation(par, source_type, laserlambda)
% This function calls the appropriate function (based on the chosen optical model) to calculate
% generation profiles as a function of position
% SOURCE_TYPE = either 'AM15' or 'laser'
% LASERLAMBDA = Laser wavelength - ignored if SOURCE_TYPE = AM15

xsolver = par.x_ihalf;
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