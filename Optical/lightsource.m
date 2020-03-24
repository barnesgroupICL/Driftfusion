function AM15 = lightsource(source_type, lambda)
% loads light sources for use with BeerLambert function
% lambda - range of wavelengths in nanometers
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
if source_type == 'AM15'
    % Load AM1.5
    AM15_data=xlsread('AM15.xls');
    AM15=1e-3*interp1(AM15_data(:,1), AM15_data(:,2), lambda, 'linear', 'extrap');
end

end