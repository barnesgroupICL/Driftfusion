function [sol_Voc, Voc] = findVocDirect(sol_ini, light_intensity, mobseti)
% Obtain approximate open circuit voltage directly using high Rs
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
sol_Voc = lightonRs(sol_ini, light_intensity, -1, mobseti, 1e6, 400);

Voct = dfana.calcVQFL(sol_Voc);
Voc = Voct(end);

end