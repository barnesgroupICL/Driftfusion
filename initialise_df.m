function initialise_df
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
% Plotting defaults
set(0,'DefaultLineLinewidth',1.5);
set(0,'DefaultAxesFontSize',16);
vert = 280;
horz = 360;
set(0,'DefaultFigurePosition', [600, 400, horz, vert]);   % [left bottom width height]
set(0,'DefaultAxesXcolor', [0, 0, 0]);
set(0,'DefaultAxesYcolor', [0, 0, 0]);
set(0,'DefaultAxesZcolor', [0, 0, 0]);
set(0,'DefaultTextColor', [0, 0, 0]);

% Add file paths to functions
addpath(genpath(pwd))


end