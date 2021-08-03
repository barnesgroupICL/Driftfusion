% An example script to demonstrate how to run a parameter exploration using
% the parallel computing toolbox
% Users will likely need to modify explore.explore2par
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
% Obtain the base parameters
par_explore = pc('input_files/3_layer_test.csv');

% For the first example we will run JV and steady-state Voc for 3 different
% active layer thicknesses and light intensities 
parex_dactive_light = explore.explore2par(par_explore, {'d(1,3)','Int'},...
    {[100e-7, 200e-7, 400e-7, 800e-7], logspace(-2,0,3)}, 200);

%% Example plots
% 1 sun JV plot for different active layer thickness
explore.plotJV(parex_dactive_light, [1,1,1,1], [0,0,1])
legend('d_{active} = 100 nm', 'd_{active} = 200 nm', 'd_{active} = 400 nm', 'd_{active} = 800 nm')

% Plot the reverse scan Voc as a function of thickness and light intensity as a surface
explore.plotsurf(parex_dactive_light, 'Voc_r', 1, 0, 0)

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder size regardless of whether it contains data or not.
% save('/Users/Username/Data/temp.mat')
